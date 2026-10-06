# Frozen package build provenance

The definitive long-format Level-1 master is generated from the frozen branch without invoking InterfaceScout prediction code.

- Branch: `benchmark-audit-freeze-2026-10-06`
- Successful GitHub Actions build run: `37537633076`
- GitHub Actions artifact: `interfacescout-benchmark-v1-2026-10-06`
- Artifact ID: `11447150593`
- GitHub artifact digest: `sha256:9bee2a265d64118a19848c69394b3c3a8c2b26bf2e543328ac580c64c16ccdeb`
- Full master rows: 910,880
- Full master columns: 44
- Duplicate source keys detected by build: 0
- Restored/verified PC-DB matrix Git blob SHA: `e04a4535579d8f092e3c61bd5681982aed20fcea`
- Payne NP metadata workbook Git blob SHA: `17167ee4cdcc77d82cf3369206a06adb049766ca`

The first package attempt exposed an integrity problem in the vendored PC-DB matrix: the large source file had been committed as an empty blob during initial text-tool vendoring. The source matrix was then restored directly from its original Chou-repository Git blob, and the destination content SHA was verified to equal the source SHA before the successful full build. No InterfaceScout output was consulted during this correction.

The user-delivered final ZIP additionally contains `benchmark_audit.xlsx` and has local SHA-256:

`532b5407f5bdf1f697d42210bb61b598df310f800b00b60433d08a9a69ec55b2`

This local delivery ZIP is a packaging layer around the frozen data artifact; benchmark eligibility and labels are unchanged.
