# InterfaceScout versioning

InterfaceScout now follows semantic versioning for public releases.

## InterfaceScout 1.0.0

The first frozen canonical release. It contains the original lightweight protein-centered compatibility model:

- side-chain relative solvent accessibility (scRSA),
- chemistry/state-conditioned local compatibility,
- C-alpha-based 5/8 A multiscale persistence,
- optional APBS potential as an auxiliary descriptor.

The immutable reference branch is `release/v1.0.0`.

## InterfaceScout 2.0.0

Current development line. This is a major version because the physical interpretation is being broadened rather than merely patched.

Planned/implemented 2.0 components include:

- biological-assembly-aware structural context,
- auxiliary CX protrusion descriptors,
- material mechanism profiles,
- applicability/limitations reporting,
- surface-type-specific physics layers,
- established nonpolar descriptors and orientation physics under development.

The current branch is `development/v2.0.0` and should be reported as `2.0.0-dev` until the physics layer and new held-out validation are completed.

## Versioning rule

- PATCH (`2.0.1`) = bug fix with no scientific-model change.
- MINOR (`2.1.0`) = backward-compatible new option/descriptor that does not change the primary model interpretation.
- MAJOR (`3.0.0`) = materially different scientific model, scoring interpretation, or incompatible API/output behavior.

Historical internal names such as `v5.1`, `v5.2`, and `v5.3` are development artifacts and are not public InterfaceScout release numbers.
