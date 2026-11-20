# Baryon YLM Superb elementals

Measurement: `BARYON_MATELEM_COLORVEC_YLM_SUPERB`.
Uses the existing Superb baryon gauge, colorvec, momentum, phase, time-selection
and batching parameters. `displacement_length` must be 1; `use_derivP` and
`displacement_list` are rejected. All operations are right derivatives.

`ylm_list` accepts explicit physical-integer descriptors only:

| Entry | Meaning |
|---|---|
| `0` | no derivatives |
| `3 m` | one derivative on quark 3, m=-1,0,1 |
| `33 L M` | two ordered derivatives on quark 3, coupled to L,M |
| `23 L M` | one derivative on quark 2 and one on quark 3, coupled to L,M |

L=0,1,2 and -L<=M<=L. No shorthand order expansion is accepted. Duplicate
entries are removed. The full basis has 22 components; all M are retained.
`drop_negative_m=true` is rejected. No vector-index permutation packing is done.

The public CG API takes **twice every spin and magnetic projection**:
`Hadron::clebsch(2,2*a,2,2*b,2*L,2*M)`. XML derivative labels use physical
integers as above. No quark spin, flavor, S3 coupling, helicity/subduction or
creation/annihilation coefficient is included in the elemental.

The Redstar baryon circular coefficients are i/2 (Nx+iNy), -i/sqrt(2) Nz,
and -i/2 (Nx-iNy) for m=+1,0,-1. For 33, CG couples d_a d_b, so Cartesian
application paths store b before a. For 23, a acts on line 2 and b on line 3.
The epsilon-contraction orientation is inherited unchanged from the existing
Superb baryon backend. Match against that Cartesian writer before comparing to
any independent convention for the Levi-Civita tensor.

Coefficients live in `baryon_derivative_ylm.cc`, not inline template headers.
The Superb prefix tree shares coordinate-space derivative intermediates across
whole vector-column blocks. Each contracted Cartesian block is converted to
double precision once, then accumulated into its coupled components. S3T storage defaults to complex float, with `storage_precision=64` available;
CG accumulation remains double and contraction precision follows the Chroma build.

## Output

S3T uses `tensorOrder=ijktdmh`, `ylm_components` for the resolved d-axis, and
`Params/ylm_list` for requests. Metadata basis:
`baryon_cg_covariant_derivatives_v1`; encoding version 1; value type 12.
The raw index tensors contain no flavor/spin permutation symmetrization.
The existing common per-quark phasing mechanism is preserved.

FileDB retains binary fields `(t_slice,left,middle,right,mom)`; left and middle
are empty, and right holds the full coupled descriptor, including placement.
Readers must explicitly support this encoding and value type. Existing FileDB
output files are rejected. S3T follows the existing writer's file creation behavior.

## Matched runtime inputs

From `tests/chroma/hadron/distillation`, use `baryon.deriv.ini.xml` and
`baryon.ylm.ini.xml` with the same gauge and colorvec inputs. They request all
22 Cartesian/coupled components respectively at 0, +/-(1,0,0), +/-(1,1,0),
+/-(1,1,1), ten vectors and four times. They produce separate `.nonzero.sdb`
files. These are S3T, not FileDB. The split Cartesian input contains ALL nine
ordered direction pairs so no permutation reconstruction is needed for comparison.

## Standalone validation

From the Chroma root, with a C++17 compiler and Python/numpy:

```sh
c++ -std=c++17 -Wall -Wextra -pedantic -Ilib \
  tests/baryon_ylm/dump_basis.cc \
  lib/meas/inline/hadron/baryon_derivative_ylm.cc \
  lib/util/ferm/clebsch.cc -o /tmp/dump_baryon_ylm
/tmp/dump_baryon_ylm > /tmp/baryon_ylm_basis.txt
python3 tests/baryon_ylm/check_basis.py /tmp/baryon_ylm_basis.txt
```

Checks cover half-integer CGs, invalid requests, deduplication, all 22 coefficient
expansions against an independent CG table, Gram normalization, and direct
circular versus Cartesian construction on a periodic noncommuting SU(3) lattice.
Finite-momentum and equal nonzero phasing cases exercise spectator antisymmetry
and the split-line exchange sign (-1)^(L+1). These are synthetic algebra tests;
actual Chroma output/readback and Redstar baryon vertex integration remain pending.

## Storage precision (27 September)

S3T now defaults to 32-bit payloads; request `storage_precision=64` for references.
FileDB remains 64-bit. The existing runtime inputs explicitly retain 64; new
`.single.ini.xml` variants select 32 and distinct output names. See
[storage validation](../ylm_storage/README.md).
