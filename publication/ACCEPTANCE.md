# Publication acceptance gates (v1.1.0-rc1)

**Historical v1.0.0 remains unchanged. Never overwrite its tag or main branch as part of candidate verification.**

## Release-blocking criteria

1. Code is on an isolated, identified source commit; Python and JavaScript dependencies and scientific model specification are recorded.
2. All cross-platform Python 3.10–3.12 unittest jobs complete successfully; no skipped numerical test is represented as passed.
3. All three OS first-install, Chromium browser, large-protein and GAPDH-tetramer jobs pass. The browser export must contain actual 3Dmol screenshots rather than fallback/generated representations.
4. Full 21-case / 231-map numerical oracle check is PASS, with all input/specification hashes valid. No acceptance thresholds are loosened in response to failures.
5. B-factor output has the expected first-model-only and selected-chain scope; a multi-model input must not copy scores onto unanalysed frames.
6. The renumbering/relabeling and insertion-code regression tests preserve all 11 chemistry values and corresponding ranked spatial regions; comparisons use identity-aligned positions.
7. Errors on invalid input and explicit warnings about excluded cofactors/glycans, missing atoms and alternate conformations are checked.
8. The final publication ZIP contains the exact source snapshot, 21 numerical JSON outputs and verification manifest, all 3 OS browser evidence ZIPs and checksum records, version notes, protocol and known limitations. The builder must reject missing or failed evidence.
9. The manuscript cites **this tested source commit** and distinguishes computed numerical verification from historical experimental comparison.

## Acceptance status

**Not yet accepted** until GitHub Actions and the strict bundle builder have passed on the **same commit**. An empty field, queued run or manually edited PASS status does not satisfy the gates.

## Distinctions

- **Label invariance:** identical coordinates/atoms after renaming or renumbering should give the same numerical surface descriptors and physically matching regions.
- **Real structural change:** a different conformation, different oligomeric assembly or omitted cofactor can genuinely change accessibility; do not force identical outputs.
- **Scientific evidence:** the frozen adsorption literature benchmark and the new code-verification panel answer different questions.
- **Coverage limit:** mmCIF-only and nonstandard complexes cannot be claimed supported simply because protein PDB tests pass.
