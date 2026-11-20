# CORR_SUPERB stored-tensor validation

The devel backend now gathers baryon tensor/derivative slots consistently, exposes the legacy averaging and zero-unsmeared settings explicitly, and retrieves reversed GenProps with gamma-five hermiticity. Insertion support uses the shorter unoriented time interval; equal source/sink endpoints retain the full diagonal GenProp tensor. Full-lattice averaging applies to products of one-point pieces, preserving mixed diagrams and exact symmetry-zero keys.

Validation on 4³×16 uses the original WEAK_FIELD/RNG, MOD64 vectors, native-double tensors, RsdCG=1e-14, and stored propagator/GenProp files. All 29 active cases passed NPT/CORR comparisons on devel and experimental no_1pt (58 comparisons, 8,380 finite values); 5 cases per chain are explicitly deferred. All Chroma jobs performed zero CG solves. rho_rho_3pt passed all 189 values per chain, maximum difference 6.3e-16. Exact keys/shapes and the unchanged per-value 1e-9 relative gate plus 1e-12 globally scaled floor were checked. Glue/tetra CORR and disconnected-loop nuc_nuc_3pt support remain deferred.

The original 12³×16 Nπ test also passed both chains: ten vectors, seven samples, 81+729 keys and 5,670 finite values. Maximum strict experimental/devel difference was 1.110668157082426e-16.

**Build qualification:** numerical executables were incremental corrected pre-YLM snapshots (baseline Chroma ad9795af plus corrected ADAT 34cddeac), not clean builds of the latest YLM checkout. The source fix is carried on devel without merging or changing the experimental graph format. The full on-demand GenProp path has separate unresolved diagnostic failures; these stored-tensor checks do not validate that path. The earlier on-demand pion two-point check passed.

See [Kanoni regression coverage](https://github.com/JeffersonLab/redstar_kanoni/blob/devel/docs/regressions.md) and the [Nπ report](https://github.com/JeffersonLab/redstar/blob/devel/docs/roper/chroma/ORIGINAL12_RUN.md). Old datasets and historical failed evidence are retained; no tolerances were weakened. no_1pt remains experimental.
