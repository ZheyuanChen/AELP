"""Validation tests for the lasy -> EPOCH raw-binary export.

These tests exercise ``Laser.write_to_file(file_format="epoch")`` from the
lasy fork (https://github.com/ZheyuanChen/lasy), which writes the complex
envelope as headerless amplitude/phase ``.dat`` files for the EPOCH
``use_spatiotemporal_profile`` injection path (epoch_dev fork).

Conventions under test (see epoch_helper.py and, on the EPOCH side,
custom_laser.f90 / laser.f90):

  lasy : E_phys = Re[E_env * exp(-i w0 t)] = |E_env| * cos(w0 t - phi)
  EPOCH: E_phys = amp * profile * SIN(integral(w dt) + phase),
         phase file consumed as-is in radians, bilinearly interpolated

  =>  profile = |E_env| / max|E_env|          (normalised to [0, 1])
      phase   = -(phi - phi_ref) + pi/2       (cos -> sin carrier shift)

  On-disk layout: C-order (n_t, n_y, n_x) doubles, i.e. Fortran
  column-major (n_x, n_y, n_t) with x fastest-varying, matching the
  epoch3d reader's ALLOCATE(matrix(n_tr1, n_tr2, n_t)).

phi_ref (envelope phase at the amplitude peak) is subtracted so that the
arbitrary global phase added by lasy's angular-spectrum propagator does not
leak into the deck; absolute carrier-envelope phase is then set by the
deck's own phase parameter. Without this referencing, two exports of the
same pulse can disagree pointwise by up to 2|sin(delta/2)| (observed as a
137% field error before the June 2025 fix).

Run with:  pytest test_lasy_epoch_export.py -v
"""

import os

import numpy as np
import pytest
from scipy.constants import c

from lasy.laser import Laser
from lasy.profiles.gaussian_profile import GaussianProfile

# Small grid: enough points to sample the carrier and the focal curvature,
# small enough that the whole suite runs in seconds.
WAVELENGTH = 0.8e-6
N_X, N_Y, N_T = 32, 24, 48
LO = (-5e-6, -4e-6, -30e-15)
HI = (5e-6, 4e-6, 30e-15)
FILE_PREFIX = "test_pulse"


@pytest.fixture(scope="module")
def exported(tmp_path_factory):
    """Build a propagated Gaussian, export it, and hand back everything the
    tests need: the lasy grid, the raw arrays reloaded from disk, and the
    metadata sidecar."""
    profile = GaussianProfile(
        wavelength=WAVELENGTH,
        pol=(1, 0),
        laser_energy=1.0,
        w0=2.0e-6,
        tau=10e-15,
        t_peak=0.0,
    )
    laser = Laser(dim="xyt", lo=LO, hi=HI, npoints=(N_X, N_Y, N_T),
                  profile=profile)
    # Propagate so the envelope picks up the angular-spectrum propagator's
    # global phase piston and a spatially varying (curved) phase front --
    # this is what makes the phi_ref referencing and the unwrapping
    # non-trivial to get right.
    laser.propagate(3e-6)

    write_dir = str(tmp_path_factory.mktemp("epoch_export"))
    laser.write_to_file(file_prefix=FILE_PREFIX, file_format="epoch",
                        write_dir=write_dir)

    shape_disk = (N_T, N_Y, N_X)
    amp = np.fromfile(os.path.join(write_dir, f"{FILE_PREFIX}_amplitude.dat"),
                      dtype=np.float64).reshape(shape_disk)
    phase = np.fromfile(os.path.join(write_dir, f"{FILE_PREFIX}_phase.dat"),
                        dtype=np.float64).reshape(shape_disk)
    meta_path = os.path.join(write_dir, f"{FILE_PREFIX}_metadata.txt")
    metadata = {}
    if os.path.exists(meta_path):
        with open(meta_path) as f:
            for line in f:
                if "=" in line:
                    key, _, rest = line.partition("=")
                    metadata[key.strip()] = rest.split()[0]
    return {
        "grid": laser.grid,
        "envelope": np.asarray(laser.grid.get_temporal_field()),  # (x, y, t)
        "amp": amp,                                   # disk order (t, y, x)
        "phase": phase,                               # disk order (t, y, x)
        "metadata": metadata,
        "write_dir": write_dir,
    }


def test_file_sizes_match_headerless_convention(exported):
    """EPOCH's only load-time sanity check is file size == product of the
    deck-declared point counts * 8 bytes; both files must satisfy it."""
    expected_bytes = N_X * N_Y * N_T * 8
    for stem in ("amplitude", "phase"):
        path = os.path.join(exported["write_dir"], f"{FILE_PREFIX}_{stem}.dat")
        assert os.path.getsize(path) == expected_bytes


def test_amplitude_normalised_to_unit_peak(exported):
    """EPOCH multiplies the file profile by the deck's amp, so the file must
    be normalised to [0, 1] with the peak exactly 1."""
    amp = exported["amp"]
    assert np.max(amp) == pytest.approx(1.0, abs=1e-15)
    assert np.min(amp) >= 0.0


def test_fortran_layout_x_fastest(exported):
    """The C-order (t, y, x) file must be the transpose of lasy's (x, y, t)
    envelope magnitude -- i.e. Fortran column-major with x fastest-varying,
    exactly matching epoch3d's ALLOCATE(matrix(n_tr1, n_tr2, n_t))."""
    env_abs = np.abs(exported["envelope"])
    expected = np.transpose(env_abs / env_abs.max(), (2, 1, 0))
    np.testing.assert_allclose(exported["amp"], expected, rtol=0, atol=1e-15)


def test_field_reconstruction_matches_lasy(exported):
    """The physics test: rebuilding E(t) the way EPOCH does,
    amp * sin(w0 t + phase), must reproduce lasy's physical field
    Re[E_env exp(-i w0 t)] (up to the deliberate phi_ref carrier
    referencing) at every grid point.

    This catches both possible sign errors: a wrong phase sign gives a
    time-reversed carrier, and a -pi/2 instead of +pi/2 carrier offset
    gives a globally sign-flipped field (a pi CEP error) -- invisible in
    intensity, but a ~200% pointwise field error against a reference run.
    """
    env = exported["envelope"]                       # (x, y, t)
    grid = exported["grid"]
    t = np.asarray(grid.axes[2])
    w0 = 2 * np.pi * c / WAVELENGTH

    # phi_ref exactly as the exporter defines it: envelope phase at the
    # amplitude peak (unwrapping only changes phi by multiples of 2*pi,
    # which cancels inside cos/sin, so wrapped angle() is fine here).
    peak_idx = np.unravel_index(np.argmax(np.abs(env)), env.shape)
    phi_ref = np.angle(env[peak_idx])

    # EPOCH-style reconstruction from the files (disk order t, y, x).
    reconstructed = exported["amp"] * np.sin(
        w0 * t[:, None, None] + exported["phase"])

    # lasy physical field with the phi_ref carrier rotation applied,
    # normalised the same way: Re[E_env e^{-i phi_ref} e^{-i w0 t}] / max.
    expected = np.real(
        env * np.exp(-1j * phi_ref) * np.exp(-1j * w0 * t)[None, None, :]
    ) / np.abs(env).max()
    expected = np.transpose(expected, (2, 1, 0))

    np.testing.assert_allclose(reconstructed, expected, rtol=0, atol=1e-9)


def test_phase_at_peak_is_carrier_offset(exported):
    """At the amplitude peak, phi - phi_ref = 0 by construction, so the
    stored phase must be exactly the +pi/2 cos->sin carrier offset."""
    peak_idx = np.unravel_index(np.argmax(exported["amp"]),
                                exported["amp"].shape)
    assert exported["phase"][peak_idx] == pytest.approx(np.pi / 2, abs=1e-12)


def test_phase_is_unwrapped_where_it_matters(exported):
    """EPOCH bilinearly interpolates the phase file, so the stored phase
    must have no 2*pi wrap seams between neighbouring samples anywhere the
    amplitude is non-negligible (interpolating across a seam would inject
    an O(pi) phase error into the boundary field)."""
    amp, phase = exported["amp"], exported["phase"]
    significant = amp > 1e-3
    for axis in range(3):
        jump = np.abs(np.diff(phase, axis=axis))
        both = np.logical_and(
            np.take(significant, range(0, significant.shape[axis] - 1),
                    axis=axis),
            np.take(significant, range(1, significant.shape[axis]),
                    axis=axis),
        )
        assert np.all(jump[both] < np.pi), (
            f"phase wrap seam along axis {axis} in a region of significant "
            "amplitude")


def test_metadata_sidecar_reports_deck_parameters(exported):
    """The binaries are headerless, so the exporter must report everything
    the deck has to declare: peak field (for amp), point counts, transverse
    extents and the time window."""
    meta = exported["metadata"]
    env_peak = np.abs(exported["envelope"]).max()
    assert float(meta["peak_field_V_per_m"]) == pytest.approx(env_peak,
                                                              rel=1e-12)
    assert int(meta["n_tr1_points"]) == N_X
    assert int(meta["n_tr2_points"]) == N_Y
    assert int(meta["n_t_points"]) == N_T
    assert float(meta["profile_tr1_min"]) == pytest.approx(LO[0])
    assert float(meta["profile_tr1_max"]) == pytest.approx(HI[0])
    assert float(meta["profile_tr2_min"]) == pytest.approx(LO[1])
    assert float(meta["profile_tr2_max"]) == pytest.approx(HI[1])
    assert float(meta["t_start"]) == pytest.approx(LO[2])
    assert float(meta["t_end"]) == pytest.approx(HI[2])


def test_rt_geometry_is_rejected(tmp_path):
    """Cylindrical grids cannot be written for EPOCH's Cartesian reader;
    the exporter must refuse loudly rather than write a misshapen file."""
    profile = GaussianProfile(
        wavelength=WAVELENGTH,
        pol=(1, 0),
        laser_energy=1.0,
        w0=2.0e-6,
        tau=10e-15,
        t_peak=0.0,
    )
    laser = Laser(dim="rt", lo=(0.0, -30e-15), hi=(5e-6, 30e-15),
                  npoints=(24, 48), profile=profile)
    with pytest.raises(NotImplementedError):
        laser.write_to_file(file_prefix="rt_pulse", file_format="epoch",
                            write_dir=str(tmp_path))
