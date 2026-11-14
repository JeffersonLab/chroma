# Apple Metal Clover SystemSolvers

Build MetalClover's public C ABI with Apple Clang, then enable the two Chroma
factories with its installation prefix:

```sh
make -C /path/to/MetalClover/metal-bsr-prototype install PREFIX=/path/to/metal-prefix
/path/to/chroma/configure --with-metal=/path/to/metal-prefix \
  --with-qdp=/path/to/qdpxx --with-datalib=/path/to/datalib [other options]
make
```

The default remains disabled (`--without-metal`). Configure checks the public
header, ABI version and library link. Chroma needs CPU FP64, Nc=3, Ns=4, four
dimensions, one QDP node, even extents and the standard QDP gamma basis.
The backend C ABI avoids passing C++ standard-library types between GCC-built
Chroma and Apple Clang's Objective-C++/libc++ Metal implementation.
Metal shader sources are embedded in the installed library.

Use `MetalBiCGClover` for Metal BiCGStab or `MetalMGClover` for flexible GMRES
with the near-null/SVD multigrid preconditioner:

```xml
<InvertParam>
  <invType>MetalMGClover</invType>
  <RsdTarget>1e-12</RsdTarget>
  <MaxIter>65536</MaxIter>
  <RHSBatch>16</RHSBatch>
  <CloverParams>
    <Mass>-0.0860</Mass>
    <u0>1</u0>
    <clovCoeffR>1.58932722549812</clovCoeffR>
    <clovCoeffT>0.902783591772098</clovCoeffT>
    <AnisoParam>
      <anisoP>true</anisoP><t_dir>3</t_dir>
      <xi_0>4.3</xi_0><nu>1.265</nu>
    </AnisoParam>
  </CloverParams>
  <Blocks>4 4 4 4</Blocks>
  <Restart>8</Restart>
  <NullVectors>24</NullVectors>
  <NullIterations>128</NullIterations>
  <NullRestart>8</NullRestart>
  <NullWarm>32</NullWarm>
  <NullBatch>4</NullBatch>
  <NullTolerance>0.001</NullTolerance>
  <SecondFilter>64</SecondFilter>
  <BottomIterations>8</BottomIterations>
  <Workers>8</Workers>
</InvertParam>
```

These are the 32-cubed header coefficients, not universal action defaults.
Use the coefficients and mass of the gauge/propagator under study. Spatial
extents divided by `Blocks`, and the temporal extent divided by its block,
must be even. The finest and second-level bases currently retain eight
vectors per chirality. The smoother is one local clover step before and
after the coarse correction; the bottom solve is sparse iterative.
Set `VECLIB_MAXIMUM_THREADS=1` to prevent nested LAPACK parallelism.
MG setup happens once in the SystemSolver constructor and is reused across
its calls, including the explicit array interface used by `SB::doInversion`: vectors of
`shared_ptr<LatticeFermion>` solutions and `shared_ptr<const LatticeFermion>`
sources. It returns one result per RHS, packs multiple RHS into shared GPU
dispatches and exposes `lastArraySolveInfo()` plus `lastBatchSolveInfo()` for
per-column and aggregate work/timing. Empty arrays are
accepted; lengths, null fields and aliased output/source fields are checked.

The adapter consumes `FermState::getLinks()`: stout smearing and fermion BC
are already present. It does not smear those links again or apply a second
temporal boundary sign. Anisotropy is folded into the backend hopping links
and into its clover construction. A CPU FP64 operator probe checks the
supplied Chroma PLUS operator against the backend before the first RHS.
The action currently requires `u0=1`; anisotropy requires `t_dir=3`; twisted mass and MaxNorm options are
rejected. The local LDL implementation requires positive local pivots and
reports failure if the action violates that assumption.

Full operators receive full RHS reduction/reconstruction; the asymmetric
odd operator receives an already reduced odd RHS for
`A_oo - D_oe A_ee^-1 D_eo`. Other subsets and symmetric-preconditioned
operators are unsupported (the operator-match probe rejects a mismatched
operator). Nonzero initial guesses and values outside the odd subset are
preserved. Each refinement residual uses the original Chroma CPU FP64
operator, and a final native check enforces `RsdTarget` on that operator.
`SystemSolverResults_t::resid` is the absolute residual, following Chroma's
contract. Relative residual and total time are printed separately.
`LinOpSysSolverMetalClover::lastSolveInfo()` exposes the detailed backend
result. The backend itself returns a result structure to its caller.

The printed rate is **nominal dispatched operator GFLOP/s**, not measured
hardware FLOP/s or total algorithm FLOP/s. The numerator includes operator,
transfer and adjoining local-vector kernels, including padded coarse
iterations. It excludes Krylov vector/reduction operations and CPU FP64
operations. The denominator includes field conversion, refinements, full
CPU FP64 checks and the final Chroma check; setup is separate. See
MetalClover's `OPERATOR_WORK.md` for formulas and scope.

The initial Chroma adapter uses in-memory FP64 spinor workspaces, alongside
Chroma's own gauge and clover allocations. The separate file-backed
standalone benchmark fitting 32³×256 does not guarantee that the full Chroma
process fits a 32 GiB machine. MPI decomposition and concurrent solves on one
context are not implemented.

## Real-gauge harness

```sh
make -C mainprogs/tests t_metal_clover
VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 \
  mainprogs/tests/t_metal_clover gauge.lime action.xml 3,3,3,8
```

The XML uses the existing perambulator metadata paths
`/DBMetaData/lattSize` and `/DBMetaData/Params/Propagator/FermionAction`.
The harness reads the real QIO/LIME gauge (including checksum checking),
reunitarizes as in gaugeStartup and constructs the native Chroma fermion
state from its XML. It tests both factory names on full and odd operators,
two RHS with one hierarchy per MG object, nonzero and converged guesses,
zero RHS, BiCG/MG solution agreement and native full reconstruction from the
odd solve. Every solve and reconstructed solution must meet 1e-12 using
original Chroma CPU FP64 operators. It is an explicit local test requiring a
gauge file, not an automatic CI test.

The regression also calls `EvenOddPrecCloverFermAct::qprop` with five full,
time-localized spin-3 sources through that exact array signature. This
exercises Chroma's vector even/odd reduction and reconstruction around the
Metal odd solver and checks every reconstructed field with the full native
FP64 operator. The distillation helpers now bound their temporary full-lattice
spinor pools by the physical spinor RHS tile and use `Ns` for the source-spin
matrix, allowing all four spin components. For a large perambulator job,
set `max_rhs` to the available host memory as well as GPU memory. The
measurement's full returned tensor still occupies caller memory; its existing
callback form can contract or write solved tiles without retaining all full
solutions. Array support does not imply that 96 full FP64 spinors fit RAM.
When SuperbBLAS is configured, the harness additionally runs the actual
`SB::ChimeraSolver` / `SB::doInversion` colorvector path used by
`inline_prop_and_matelem_distillation_superb_w.cc`: five deterministic complex
colorvectors, spin 3, one time source and `max_rhs=3` (GPU tiles of 3+2). The extracted full
solutions must pass the original Chroma FP64 operator at 1e-12. This checks
the tensor-to-source spin convention as well as the array SystemSolver call.

## Initial integration validation

Before the batched GPU extension, a fresh FP64 Chroma build on the requested `rge/devel-datalib-reorg` base
passed both Metal-enabled and Metal-disabled configure checks. The complete
enabled library and programs built with MPI QDP++/SuperbBLAS on one rank.
The real 12³×96 LIME field at mass −0.0743 passed all 20 scalar, native
reconstruction, qprop-array and actual SuperbBLAS distillation checks. The
maximum original Chroma CPU FP64 relative residual was 2.3534e-14 and the
maximum whole-lattice backend/native operator probe difference was
1.8249e-16. Each MG object reused one hierarchy for its two sources. These
are correctness checks, not large-volume Chroma performance measurements.
The structured results and input/source hashes are recorded in MetalClover's
`results/near-null-svd32/chroma-integration/` evidence directory.

The SuperbBLAS build also supports disabling PRIMME. Its optional Generalized
Davidson helper reports that `--with-primme` is required when that helper is
requested; Metal multigrid does not require PRIMME.

## Batched GPU RHS

`RHSBatch` is the maximum resident GPU RHS (1..64, default 8). The library
may reduce it to respect its estimated device-memory budget. Cap 1 uses the
scalar path. It is independent of MPI geometry, CPU `Workers`, and the
measurement's `max_rhs`: `max_rhs` controls the array delivered to qprop and
bounds its native CPU temporary fields, while `RHSBatch` controls device tiles.
Gauge, clover, transfer and coarse-matrix data are shared by all widths; changing
width allocates RHS workspaces without rebuilding the hierarchy. Fine/coarse
products, smoothers, transfers and bottom iterations all run with the RHS
width. Krylov coefficients, convergence and FP64 acceptance remain independent
for each column, including zero, collinear and initially converged columns.

Each correction residual is normalized independently before FP32 conversion
and its correction is rescaled in FP64. All outer native Schur residuals and
final native full-qprop residual evaluations remain CPU FP64. Per-column solve
times are amortized shares of the batch time; detailed aggregate batch work,
actual maximum width, tile count and workspace time are exposed separately.
The printed work model includes inactive/padded lanes and excludes CPU FP64
and Krylov-vector arithmetic from its numerator, as documented in MetalClover.

The real-field throughput harness compares RHS tiles 1, 8 and 16 on the same
sources using one qprop object/hierarchy per method throughout, then doubles
the source collection. It checks every full original Chroma FP64 residual at
1e-12 and compares the first collection's complete solutions against a scalar
reference. Source/solution storage and reference-file I/O are outside solve
timing; native qprop reduction, reconstruction and full residual evaluation
are inside. An additional independent full check is performed after timing.

```sh
make -C mainprogs/tests t_metal_clover_batch
VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 \
  mainprogs/tests/t_metal_clover_batch gauge.lime action.xml /tmp/batch12 24
```

The current harness takes BASE_NUM_VECS: 24 vectors mean 96 physical RHS
in each of the three modes, followed by 48 vectors / 192 RHS in tiles of 16.
It is a single-rank test for the real 12³×96 action XML, with one scalar
reference file and CSV output. Choose a fresh output prefix. The build requires
MetalClover's additive batch API (`METALCLOVER_BATCH_API_VERSION >= 1`);
configure verifies the header and both new batch-control/info symbols.

The batched regression subsequently passed 32 native scalar/reconstruction/
qprop/Superb checks at 1e-12, with maximum measured full relative residual
2.5474e-14. Both inverters used actual GPU tiles of three plus two for the
five-source array tests. Separate real-field C ABI tests also cover mixed
source scales, zero/exact/nonzero guesses, collinear RHS, partial tiles,
iteration budgets and scalar agreement; a GCC-built client passes against the
Apple Clang backend. MetalClover's `BATCH_RHS.md` records the 96–192 RHS native
qprop throughput measurements and the detailed counting/timing scope.

The original batch timings used independent full-field random spinors;
they are physical RHS counts, rather than color eigenvectors. On the real
12³×96 mass −0.0743 field, the same 96 sources took 3.358,
2.295 and 2.201 s/RHS for BiCGStab at GPU widths 1, 8 and 16; MG took
2.488, 1.416 and 1.382 s/RHS. Width sixteen improved throughput by 1.53×
and 1.80× respectively without changing iteration counts. The additional
192-source runs took 2.214 and 1.395 s/RHS at width sixteen. These are
single wall-clock measurements of complete native qprop calls, including
FP64 checks and excluding setup. The one reused MG hierarchy took 19.714 s
to construct. All full native residuals were below 2.433e-14, and the first
96 full solutions agreed across methods/widths within 6.343e-13 relative.
Peak process RSS was 14.547 GiB on the 32 GiB M1 Max (24 GPU cores).

## Eloy's distillation caller with a batched Metal backend

The production measurement retains Eloy's original source-spin outer loop,
eigenvector batching within each spin, contractions and storage organization.
`num_vecs` eigenvectors still require `Ns * num_vecs` physical inversions per
source time and phasing. That total does not require coalescing all source
spins in one call. One solver object and MG hierarchy serve all the loops.

`max_rhs` is the physical spinor count in a caller batch; with a fixed source
spin it also counts the eigenvectors in that call. A positive cap is used
without rounding, and **max_rhs=1 is valid**. For three vectors at cap two,
each spin sends batches of two and one. For 5,184 vectors at cap sixteen,
four spins produce 1,296 calls and 20,736 inversions. There is no small cap
on production `num_vecs`. The Metal `RHSBatch` is a separate device/memory cap;
set both to sixteen for the largest width in the existing measured comparison.

The original perambulator measurement and header are restored byte-for-byte
from the pre-overhaul branch. Its helper retains the localized bounded
source/solution pools and `Ns` spin dimension fix. Metal's true array backend,
shared gauge/clover/hierarchy data, per-column convergence and full original
Chroma CPU FP64 residual checks remain in place.

The measurement contracts and stores one source-spin channel and sink phasing
at a time, removing the overhaul's simultaneous four-spin projected buffers.
Sink eigenvector projection can reduce retained data when the basis is
truncated. `max_rhs` bounds inversion temporaries; the returned tensor,
eigenvectors and projected output still consume memory. The existing callback
interface remains available for a separate streaming unsmeared-output writer.

The throughput harness keeps BASE_NUM_VECS and explicit physical-RHS reporting,
but now follows the production spin-outer ordering at source time zero. It
compares caller caps 1, 8 and 16 and doubles the vector count in the last mode.
A single vector therefore makes four one-RHS calls even with a larger cap.
Its 24/48-vector memory bound is test-only, not a production restriction.
The previous coalescing test evidence at commit `bc151f436` remains historical.

The restored real-gauge regression passed all 64 native checks at source time
zero. Both inverters respected `max_rhs=1` (four separate one-RHS calls for
one eigenvector) and `max_rhs=2` (2+1 batches for each spin of three vectors).
The new spin-source cases had maximum full original Chroma FP64 residual
2.0404e-14 and projected-element difference 8.882e-16 against direct CPU
inner products. The updated spin-major harness smoke test also passed its
full-solution comparisons with maximum residual 2.4445e-14 and BiCGStab/MG
difference 2.5393e-13. At that regression stage the full inline file-I/O and full
5,184-vector jobs were not run; the separate streaming job is recorded below. MetalClover's `results/eloy-spin-batching12/` contains the
structured evidence and exact source/binary hashes.


## Streaming unsmeared time-zero output

The separate `t_metal_clover_unsmeared` program uses the existing public
`SB::doInversion` callback. It leaves Eloy's projected perambulator measurement
and helper unchanged. Build it explicitly with MetalClover and SuperbBLAS enabled:

```sh
make -C mainprogs/tests t_metal_clover_unsmeared
mainprogs/tests/t_metal_clover_unsmeared \
  gauge.lime action.xml eigs.t_0.mod1000 /path/to/output-prefix \
  0 5184 16 MetalMGClover
```

Arguments after the prefix are first vector, vector count, physical caller RHS
cap and method (`MetalMGClover` or `MetalBiCGClover`). The current program targets
one rank, 12³×96, source and sink at time zero. Four source spins for 5,184
color eigenvectors require 20,736 physical inversions. The original spin-major
loop and one solver/hierarchy span the whole collection; cap one is supported.
The device cap is sixteen. Stored FP32 MOD eigenvectors are promoted to FP64
without changing ordering or phases.

Every full solution must pass both the original native qprop FP64 residual
and a separate `UnprecCloverLinOp` CPU FP64 evaluation at relative tolerance
1e-12. The latter checks all 96 sink times before retaining only time zero.
The output is point sinks, without sink eigenvector projection. The version-1
NumPy file is complex128, C order, with axes
`(source_vector, source_spin, z, y, x, sink_colour, sink_spin)` and shape
`(5184,4,12,12,12,3,4)`. Spatial order is natural x-fastest order, converted
explicitly from QDP layout. The complete file occupies about 6.41 GiB.

The `.meta.xml` sidecar includes action/source metadata and file hashes;
`.columns.csv` journals residuals, source norms, iterations and column checksums.
Payload space is reserved before inversion. A `.state` checkpoint is published
atomically after syncing every accepted batch. SIGTERM/SIGINT stops after the
current batch with exit 75. Resume with the same arguments plus `--resume`;
immutable metadata and all committed payload checksums are verified first.
An already complete resume only verifies the output and skips MG setup.
`--stop-after-batches N` supports bounded recovery tests.

MetalClover's `metal-bsr-prototype/run_unsmeared.py` adds process/disk monitoring,
idle-sleep prevention, session logs, final native checksum rereading, NPY/journal
validation and a whole-file SHA256 completion report. Its `UNSMEARED_T0.md`
documents the output reader and exact real-field production inputs. The large
output and eigenvector files are not committed.


The full streaming run on configuration 1000 completed on 2026-11-14. All
5,184 supplied t=0 eigenvectors and four source spins (20,736 physical RHS)
passed both full original Chroma CPU FP64 residual checks at 1e-12. The maximum
relative residual was 2.4063810813519627e-14. One solver/hierarchy served all
1,296 width-sixteen array calls in one uninterrupted session. The subsequent
completed-output resume verifier reread all 20,736 column checksums successfully.

The 6,879,707,264-byte complex128 point-sink NPY file is retained at
`/Users/admin/Downloads/MetalClover-results/cfg1000-t0-unsmeared.npy`.
Its SHA256 is
`3b017ee3fd6989514f3ac23f79859082442da18b3fb201c278d14fd85c8ffc51`.
Total solve-session wall time, including input, setup, both residual checks and
streamed output, was 30,140.426 s (8 h 22 min 20 s). Peak process RSS was
5.955 GiB on an M1 Max with 24 GPU cores and 32 GiB memory. Nominal dispatched
GPU operator work divided by this wall time was 277.49 GFLOP/s; CPU FP64 and
Krylov-vector arithmetic are excluded from that numerator. Final checksum
rereading and SHA256 are outside this timing. MetalClover's
`results/unsmeared-t0-5184/` records the complete column journal, metadata,
logs, resource samples, frozen hashes and production validation. The original
projected perambulator measurement and helper were unchanged for this addition.
