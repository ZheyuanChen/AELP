This cell of the design (amplitude via native deck expression, phase via
binary file) is **not reachable with the current epoch3d code** and is
intentionally left empty.

Reason (found while building this test, epoch3d `custom_laser.f90`,
`custom_laser_spatial_setup` / `load_spatial_fields`, commit `537b0445`):
on the static spatial path, `use_phase_from_file` is only ever consulted
*inside* `load_spatial_fields`, which itself only runs when
`use_custom_profile = T` -- and when it runs, it unconditionally loads
`profile_data_file` into `laser%profile` first. There is no deck
combination that gets phase from a file while leaving amplitude on the
native deck-expression path; `use_custom_profile = T` always pulls
amplitude from a file too.

This reduces the achievable design from a true 2x2 to a 3-cell chain:
  1. `amp_deck_phase_deck/` -- baseline (fully native, ground truth)
  2. `amp_file_phase_deck/` -- isolates the AMPLITUDE injector (phase held
     on the native deck expression)
  3. `amp_file_phase_file/` -- both from file (total combined effect;
     compare against cell 2 to isolate the PHASE injector's additional
     contribution)

See `../README.md` for the full reasoning and `../analyse.py` for how the
three cells are compared. Worth flagging to a future epoch_dev session as
a possible small enhancement (decoupling the two gates) if a use case
needs amplitude-native + phase-file specifically (e.g. LASY phase with a
synthetic/deck amplitude).