"""Validation tests for the lasy -> EPOCH raw-binary export.

These tests exercise ``Laser.write_to_file(file_format="epoch")`` (3D,
epoch3d spatiotemporal path) and ``file_format="epoch2d"`` (single lasy-y
plane, epoch2d spatiotemporal path) from the lasy fork
(https://github.com/ZheyuanChen/lasy), which writes the complex envelope
as headerless amplitude/phase ``.dat`` files for the EPOCH
``use_spatiotemporal_profile`` injection path (epoch_dev fork).

Conventions under test (see epoch_helper.py and, on the EPOCH side,
custom_laser.f90 / laser.f90):

  lasy : E_phys = Re[E_env * exp(-i w0 t)] = |E_env| * cos(w0 t - phi)
  EPOCH: E_phys = amp * profile * SIN(integral(w dt) + phase),
         phase file consumed as-is in radians, bilinearly interpolated

  =>  profile = |E_env| / max|E_env|          (normalised to [0, 1])
      phase   = -(phi - phi_ref) + pi/2 + carrier_phase_ref

  On-disk layout, 3D: C-order (n_t, n_y, n_x) doubles, i.e. Fortran
  column-major (n_x, n_y, n_t) with x fastest-varying, matching the
  epoch3d reader's ALLOCATE(matrix(n_tr1, n_tr2, n_t)).
  On-disk layout, 2D: C-order (n_t, n_x) doubles, i.e. Fortran
  (n_transverse, n_t) transverse fastest, matching the epoch2d reader's
  ALLOCATE(file_field_matrix(n_transverse_points, n_t_points)); the
  retained lasy x axis maps to EPOCH2d's transverse y axis.

phi_ref (envelope phase at the amplitude peak) is subtracted so that the
arbitrary global phase added by lasy's angular-spectrum propagator does not
leak into the deck; the absolute carrier-envelope phase is then set by
``carrier_phase_ref`` (default 0, i.e. stored phase at the peak is exactly
+pi/2). Without this referencing, two exports of the same pulse can
disagree pointwise by up to 2|sin(delta/2)| (observed as a 137% field
error before the June 2025 fix).

Run with:  pytest test_lasy_epoch_export.py -v
"""

import os

import numpy as np
import pytest
from scipy.constants import c

from lasy.laser import Laser
from lasy.profiles.gaussian_profile import GaussianProfile

# Small grid: enough points to sample the carrier and the focal curvature,
# small enough that the whole suite runs in seconds. N_Y is odd so the
# y = 0 plane lies exactly on the grid (used by the 2D-slice cross-checks).
WAVELENGTH = 0.8e-6
N_X, N_Y, N_T = 32, 25, 48
LO = (-5e-6, -4e-6, -30e-15)
HI = (5e-6, 4e-6, 30e-15)
FILE_PREFIX = "test_pulse"


def make_laser():
    """A propagated Gaussian: propagation adds the angular-spectrum
    propagator's global phase piston and a spatially varying (curved)
    phase front, which is what makes the phi_ref referencing and the
    unwrapping non-trivial to get right."""
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
    laser.propagate(3e-6)
    return laser


def read_export(write_dir, shape):
    """Reload an exported amplitude/phase/metadata triple from disk."""
    amp = np.fromfile(os.path.join(write_dir, f"{FILE_PREFIX}_amplitude.dat"),
                      dtype=np.float64).reshape(shape)
    phase = np.fromfile(os.path.join(write_dir, f"{FILE_PREFIX}_phase.dat"),
                        dtype=np.float64).reshape(shape)
    metadata = {}
    meta_path = os.path.join(write_dir, f"{FILE_PREFIX}_metadata.txt")
    with open(meta_path) as f:
        for line in f:
            if "=" in line and not line.lstrip().startswith("#"):
                key, _, rest = line.partition("=")
                metadata[key.strip()] = rest.split()[0]
    return amp, phase, metadata


@pytest.fixture(scope="module")
def exported(tmp_path_factory):
    """Build one laser, export it in both EPOCH formats, and hand back the
    envelope plus the raw arrays reloaded from disk."""
    laser = make_laser()
    env = np.asarray(laser.grid.get_temporal_field())  # lasy order (x, y, t)

    dir3d = str(tmp_path_factory.mktemp("epoch_export_3d"))
    laser.write_to_file(file_prefix=FILE_PREFIX, file_format="epoch",
                        write_dir=dir3d)
    amp3, phase3, meta3 = read_export(dir3d, (N_T, N_Y, N_X))

    dir2d = str(tmp_path_factory.mktemp("epoch_export_2d"))
    laser.write_to_file(file_prefix=FILE_PREFIX, file_format="epoch2d",
                        write_dir=dir2d)
    amp2, phase2, meta2 = read_export(dir2d, (N_T, N_X))

    return {
        "grid": laser.grid,
        "envelope": env,
        "amp": amp3, "phase": phase3, "metadata": meta3, "write_dir": dir3d,
        "amp2d": amp2, "phase2d": phase2, "metadata2d": meta2,
        "write_dir2d": dir2d,
    }


# ---------------------------------------------------------------------------
# 3D path (epoch3d spatiotemporal)
# ---------------------------------------------------------------------------

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
    extents and the time window (epoch3d element names)."""
    meta = exported["metadata"]
    env_peak = np.abs(exported["envelope"]).max()
    assert float(meta["peak_field_V_per_m"]) == pytest.approx(env_peak,
                                                              rel=1e-12)
    assert int(meta["n_tr1"]) == N_X
    assert int(meta["n_tr2"]) == N_Y
    assert int(meta["n_t_points"]) == N_T
    assert float(meta["tr1_min"]) == pytest.approx(LO[0])
    assert float(meta["tr1_max"]) == pytest.approx(HI[0])
    assert float(meta["tr2_min"]) == pytest.approx(LO[1])
    assert float(meta["tr2_max"]) == pytest.approx(HI[1])
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
    for fmt in ("epoch", "epoch2d"):
        with pytest.raises(NotImplementedError):
            laser.write_to_file(file_prefix="rt_pulse", file_format=fmt,
                                write_dir=str(tmp_path))


# ---------------------------------------------------------------------------
# 2D path (epoch2d spatiotemporal, lasy-y slice)
# ---------------------------------------------------------------------------

def test_2d_file_sizes_match_headerless_convention(exported):
    """epoch2d declares n_transverse_points * n_t_points; file size must be
    that times 8 bytes."""
    expected_bytes = N_X * N_T * 8
    for stem in ("amplitude", "phase"):
        path = os.path.join(exported["write_dir2d"],
                            f"{FILE_PREFIX}_{stem}.dat")
        assert os.path.getsize(path) == expected_bytes


def test_2d_layout_matches_y0_slice(exported):
    """The C-order (t, x) file must be the transpose of the y = 0 plane of
    lasy's (x, y, t) envelope magnitude, normalised to the slice's own
    peak -- Fortran (n_transverse, n_t) with the transverse axis fastest,
    matching epoch2d's ALLOCATE(file_field_matrix(n_transverse, n_t))."""
    j0 = N_Y // 2  # odd N_Y, symmetric bounds -> y = 0 exactly
    env_slice = np.abs(exported["envelope"][:, j0, :])
    expected = np.transpose(env_slice / env_slice.max(), (1, 0))
    np.testing.assert_allclose(exported["amp2d"], expected, rtol=0,
                               atol=1e-15)
    assert np.max(exported["amp2d"]) == pytest.approx(1.0, abs=1e-15)


def test_2d_field_reconstruction_matches_lasy(exported):
    """Same pointwise physics check as the 3D version, on the slice:
    amp * sin(w0 t + phase) rebuilt from the 2D files must reproduce
    lasy's physical field in the y = 0 plane (up to the phi_ref carrier
    referencing, which for the 2D export is defined on the slice)."""
    j0 = N_Y // 2
    env = exported["envelope"][:, j0, :]             # (x, t)
    t = np.asarray(exported["grid"].axes[2])
    w0 = 2 * np.pi * c / WAVELENGTH

    peak_idx = np.unravel_index(np.argmax(np.abs(env)), env.shape)
    phi_ref = np.angle(env[peak_idx])

    reconstructed = exported["amp2d"] * np.sin(
        w0 * t[:, None] + exported["phase2d"])
    expected = np.real(
        env * np.exp(-1j * phi_ref) * np.exp(-1j * w0 * t)[None, :]
    ) / np.abs(env).max()
    expected = np.transpose(expected, (1, 0))

    np.testing.assert_allclose(reconstructed, expected, rtol=0, atol=1e-9)


def test_2d_phase_at_peak_is_carrier_offset(exported):
    """Stored 2D phase at the slice's amplitude peak must be exactly +pi/2
    (default carrier_phase_ref = 0)."""
    peak_idx = np.unravel_index(np.argmax(exported["amp2d"]),
                                exported["amp2d"].shape)
    assert exported["phase2d"][peak_idx] == pytest.approx(np.pi / 2,
                                                          abs=1e-12)


def test_2d_phase_is_unwrapped_where_it_matters(exported):
    """No 2*pi wrap seams in the 2D phase file where the amplitude is
    non-negligible (epoch2d bilinearly interpolates in (y, t) too)."""
    amp, phase = exported["amp2d"], exported["phase2d"]
    significant = amp > 1e-3
    for axis in range(2):
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


def test_2d_consistent_with_3d_export(exported):
    """The 2D export must be exactly the y = 0 plane of the 3D export: for
    a beam whose global peak lies in that plane, the normalisation constant
    and phi_ref are identical, so amplitude and phase agree sample by
    sample (no separate conversion drift between the two paths)."""
    j0 = N_Y // 2
    np.testing.assert_allclose(exported["amp2d"], exported["amp"][:, j0, :],
                               rtol=0, atol=1e-15)
    # Phase may legitimately differ by exact multiples of 2*pi (the branch
    # anchoring runs over different domains); compare modulo 2*pi.
    dphi = exported["phase2d"] - exported["phase"][:, j0, :]
    dphi_mod = np.abs(dphi - 2 * np.pi * np.round(dphi / (2 * np.pi)))
    significant = exported["amp2d"] > 1e-3
    assert np.all(dphi_mod[significant] < 1e-12)


def test_2d_metadata_reports_slice_and_deck_parameters(exported):
    """The 2D sidecar must record the actual slice coordinate and the
    epoch2d deck elements (n_y from lasy's x axis, which becomes EPOCH2d's
    transverse y)."""
    meta = exported["metadata2d"]
    j0 = N_Y // 2
    env_slice_peak = np.abs(exported["envelope"][:, j0, :]).max()
    assert float(meta["peak_field_V_per_m"]) == pytest.approx(env_slice_peak,
                                                              rel=1e-12)
    assert float(meta["slice_coord_actual"]) == pytest.approx(0.0, abs=1e-20)
    assert int(meta["n_y"]) == N_X
    assert int(meta["n_t_points"]) == N_T
    assert float(meta["y_min"]) == pytest.approx(LO[0])
    assert float(meta["y_max"]) == pytest.approx(HI[0])
    assert float(meta["t_start"]) == pytest.approx(LO[2])
    assert float(meta["t_end"]) == pytest.approx(HI[2])


def test_carrier_phase_ref_shifts_stored_phase(exported, tmp_path):
    """carrier_phase_ref must add exactly that constant to the stored
    phase -- it is the only way to set the carrier-envelope reference,
    since EPOCH ignores the deck's phase expression when the phase comes
    from file (e.g. the Gouy-phase pin used when comparing against a
    native-deck run)."""
    ref = -0.5536  # a Gouy-like pin, half atan(2), as used in campaign B
    laser = make_laser()
    laser.write_to_file(file_prefix=FILE_PREFIX, file_format="epoch2d",
                        write_dir=str(tmp_path), carrier_phase_ref=ref)
    _, phase_shifted, meta = read_export(str(tmp_path), (N_T, N_X))
    np.testing.assert_allclose(phase_shifted - exported["phase2d"], ref,
                               rtol=0, atol=1e-12)
    assert float(meta["carrier_phase_ref"]) == pytest.approx(ref)
