# CORR_SUPERB one-point time behavior

On `devel-rge-no_1pt`, CORR_SUPERB requires the following element in its Param:

```xml
<average_1pt_diagrams>false</average_1pt_diagrams>
```

False evaluates ordinary version-4 graphs within their requested time windows,
including disconnected one-point pieces. It does not expand them over all lattice
time slices. Multiple sources retain their distinct graph time descriptors.

True is currently unsupported. Reading it raises
`CORR_SUPERB: average_1pt_diagrams=true is unsupported` before contraction work.
Omitting the element is an input error. Parameter output includes the flag.
Programmatically constructed measurements and the ADAT backend also reject true.

Use matching ADAT/redstar-datalib and Redstar feature branches and regenerate
older graph files. The shared graph format is 4. Redstar_npt separately supports
both averaging modes; support for true in CORR_SUPERB remains future work.

Validation: the parameter-reader harness accepts false and rejects true and missing input. The corrected original 12³×16 Nπ execution passed both chains with ten vectors, seven samples, 81+729 keys and 5,670 finite values. Strict experimental/devel maximum difference was 1.110668157082426e-16. The compact earlier evidence is retained.

The separate native-double 4³×16 RsdCG=1e-14 stored-tensor regression matrix passed all 29 active cases on each chain, with five explicit CORR deferrals per chain. All 189 rho_rho_3pt values passed (maximum difference 6.3e-16) after correcting reversed and diagonal GenProp lookup. All Chroma comparison jobs read stored perambulators/GenProps and performed zero CG solves. The source also backports devel meson callback ordering, partial propagator batches and reversed-time fallback handling.

Numerical binaries were incremental corrected pre-YLM snapshots; this does not claim a clean latest-YLM build, MPI or GPU validation. Separate full on-demand GenProp diagnostics remain unresolved and are not covered by these cached-tensor passes. no_1pt remains experimental, true remains unsupported, and format 4 is distinct from the devel regression baseline. See the [Nπ correction](https://github.com/JeffersonLab/redstar/blob/devel/docs/roper/chroma/NPI_FIX.md) and [Kanoni coverage](https://github.com/JeffersonLab/redstar_kanoni/blob/devel/docs/regressions.md).
