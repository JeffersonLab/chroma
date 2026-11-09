# Standalone datalib storage and CORR_SUPERB

This reorg branch requires `--with-datalib=PREFIX` for shared data and storage.
Add `--with-superbblas=PREFIX` to enable the SUPERB measurements, including
CORR_SUPERB. Use the same compiler/OpenMP runtime and matching QDP/QMP/SuperbBLAS
builds. No ADAT installation is required.

Legacy measurement writers now use `Datalib::TensorStorage::StorageMeson`,
`StorageBaryon`, `StorageProp4`, `StorageGenprop4`, `StorageDisco`, `StorageGlue`
or the common `StorageRecords` interface. Their record operations preserve
existing encoded keys, values and metadata. Persistence and backend selection
belong to datalib. The temporary `filedb_qdp.h` and `key_val_db.h` headers and
Chroma database/serializable wrapper classes are removed.

`util/ferm/storage_io.h` contains only a QDP local serialization helper.
Chroma owns rank selection, file-existence broadcasts and error/abort handling.
Collective producers perform the existing reductions on all ranks, then encode
and store records on the primary rank. Local SUPERB producer paths keep their
existing ownership. Datalib performs no QDP or MPI communication.

QDP XML and lattice I/O remain collective; datalib XML I/O is serial.
`SB::StorageTensor` continues to use the consumer's SuperbBLAS communicator.
Chroma obtains the storage component's compile/link flags through
`datalib-config --storage-cxxflags`, `--storage-ldflags` and `--storage-libs`.
Both Autotools and CMake use those backend-neutral query names. With tensor
support enabled, Chroma and its clients use the matching `USE_SUPERBBLAS`
profile because the storage class layouts depend on it. `chroma-config` exports
that profile and the matching SuperbBLAS includes.

See datalib's `docs/reorg/CHROMA-STORAGE-REGRESSIONS.md` for the tested
Autotools stack, old/new writer comparisons, full native-double matrix and
retained deferrals and limits. Full multi-rank CORR, GPU and on-demand GenProp
coverage are outside the cached-tensor suite.
