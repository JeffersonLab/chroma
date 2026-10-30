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

Validation: the parameter-reader harness accepts false and rejects true and
missing input. The complete translation unit passes a compilation check against
the matching feature ADAT headers. A full Chroma numerical execution has not yet
been validated.
