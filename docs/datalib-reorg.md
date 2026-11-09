# Standalone datalib and CORR_SUPERB

This reorg branch requires `--with-datalib=PREFIX` for shared data and FILEDB
support. Add `--with-superbblas=PREFIX` to enable the SUPERB measurements,
including CORR_SUPERB. Both BUILD_SB and BUILD_REDSTAR_DATALIB must be defined.
Use the same compiler/OpenMP runtime and matching QDP/QMP/SuperbBLAS builds.
No ADAT installation is required.

Datalib owns FILEDB and the typed meson/baryon/propagator/GenProp/DISCO storage
interfaces. The older measurement code is unchanged. Its former QDP FILEDB
adapter is retained as `util/ferm/filedb_qdp.h` until that code is modernized.
The data-parallel XML boundary and SuperbBLAS communicator selection are unchanged.

Datalib's `docs/reorg/FILEDB-TOOLS-REGRESSIONS.md` records the full native-double
paired suite, retained deferrals, clean source provenance and platform-test limits.
