# Single-precision YLM storage

New meson and baryon YLM Superb writers default to `<storage_precision>32</storage_precision>`
for S3T. Set 64 for double-precision references. Cartesian writers are unchanged.
Coordinate-space contraction follows the Chroma build precision; CG accumulation
remains double, with rounding only when the completed tensor is saved.

FileDB currently requires 64 (its default). Requesting 32 with
`use_superb_format=false` is rejected rather than changing the old binary value
schema silently. Metadata identifies basis, hadron family and payload precision.
The S3T binary dtype is authoritative and must agree with metadata.

Existing YLM test inputs now explicitly request 64. New `*.single.ini.xml`
variants request 32 and write distinct `.single.sdb` outputs, preserving references.
No new Chroma runtime output is included in this change.

Standalone round-trip test (requires Superb headers):

```sh
c++ -std=c++17 -I/path/to/superbblas/include \
  tests/ylm_storage/roundtrip.cc -o /tmp/ylm-storage-roundtrip
/tmp/ylm-storage-roundtrip /tmp/ylm-storage-roundtrip.s3t
```

This exercises double-to-float S3T writing, binary dtype, checksums and float-to-double
readback. It does not validate Chroma contractions, ADAT readers or ColorVec nodes.
ADAT/ColorVec YLM integration is still pending on the corresponding test branches.
