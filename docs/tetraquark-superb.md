# Tetraquark elemental SUPERB storage

`TETRAQUARK_MATELEM_COLORVEC` now defaults to native SUPERB tensor output.
The code ports the original Harom tetraquark construction into native Chroma;
the two color-coupling formulas are preserved. Set `use_superb_format` to false
inside the measurement Param only when explicitly writing legacy SDB output.
Use an explicit `.sbt` filename for tensor output. The filename is not a converter.

Schema version 1 uses complex<double> values and rank-seven `ijkltmc` axes:
four vector indices (i fastest), time, momentum and color construction.
Metadata records the tensor order, lattice, decay direction, producer precision,
parameters, explicit momentum dictionary and the two color constructions
`3barx3->1,1` and `6x6bar->1,1`. Sparse blocks have block checksums.
Native double MOD color vectors are read without a float staging buffer.
The source metadata XML declaration is removed before embedding it in Source_info,
so the Chroma output remains well-formed XML.

ADAT `io/tetraquark_superb.h` exposes a typed native reader with legacy SDB
compatibility. ColorVec's DBTetraQuarkElementalOperator now uses it, retaining
existing cache and physics interfaces. Native input requires USE_SUPERBBLAS.
The reader validates schema, payload type, axes and dictionaries, and rejects
missing or incomplete vector blocks rather than treating them as zero.
Rebuild ColorVec and then relink Redstar after installing the ADAT header.

The initial native 4³×16 check used ten double-precision vectors, all sixteen
times, both color constructions and zero momentum: 32 keys / 320,000 finite
complex values agreed exactly with legacy SDB output (maximum absolute error 0).
Missing color, momentum, out-of-range time and an intentionally incomplete vector block were rejected. The full tetra_F27_S2 Redstar regression also passed through the rebuilt
ColorVec consumer: three correlator keys / 48 finite complex values at tolerance
1e-12, with byte-identical legacy/native output databases. Both runs use the same
regenerated double-precision propagator input and unchanged graph/operators. The test executable is an incremental rebuild of the recorded pre-YLM
Chroma snapshot; it does not establish a clean latest-YLM build.

The historical SDB dataset and regression baselines remain separate. no_1pt
remains experimental; this storage change does not change corr_graph.bin.

Without BUILD_SB, the measurements retain legacy output by default and reject an
explicit request for SUPERB output. Both configurations compile.
