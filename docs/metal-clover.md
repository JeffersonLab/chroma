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
sources. It returns one result per RHS, processes one GPU RHS at a time and
exposes `lastArraySolveInfo()` for per-column work/timing. Empty arrays are
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

The regression also calls `EvenOddPrecCloverFermAct::qprop` with two full,
time-localized spin-3 sources through that exact array signature. This
exercises Chroma's vector even/odd reduction and reconstruction around the
Metal odd solver and checks every reconstructed field with the full native
FP64 operator. The distillation helpers now bound their temporary full-lattice
spinor pools by `min(num_vectors,max_rhs)` and use `Ns` for the spin-source
matrix, allowing the fourth spin component. For a large perambulator job,
set `max_rhs` to the available host memory as well as GPU memory. The
measurement's full returned tensor still occupies caller memory; its existing
callback form can contract or write solved tiles without retaining all full
solutions. Array support does not imply that 96 full FP64 spinors fit RAM.
When SuperbBLAS is configured, the harness additionally runs the actual
`SB::ChimeraSolver` / `SB::doInversion` colorvector path used by
`inline_prop_and_matelem_distillation_superb_w.cc`: two deterministic complex
colorvectors, spin 3, one time source and `max_rhs=1`. The extracted full
solutions must pass the original Chroma FP64 operator at 1e-12. This checks
the tensor-to-source spin convention as well as the array SystemSolver call.

## Local validation

A fresh FP64 Chroma build on the requested `rge/devel-datalib-reorg` base
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
