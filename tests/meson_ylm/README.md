# Coupled meson derivative elementals

New measurement: `MESON_MATELEM_COLORVEC_YLM_SUPERB`.
Use the existing superb meson measurement's gauge, colorvec, smearing,
time/momentum and output parameters. Example Param fragment:

```xml
<ylm_list>
  <elem>0</elem>
  <elem>1 -1</elem>
  <elem>1 0</elem>
  <elem>1 1</elem>
  <elem>2 0 0</elem>
  <elem>3 1 0 0</elem>
  <elem>3 2 3 0</elem>
</ylm_list>
<drop_negative_m>false</drop_negative_m>
<displacement_length>1</displacement_length>
<use_superb_format>false</use_superb_format>
```

YLM elementals always use derivatives. The obsolete `use_derivP` input tag
is rejected; remove it from existing inputs.

This example selects seven components, not the complete basis.
The complete 40-component input is `numerics/inputs/meson.ylm.ini.xml`
(relative to the notes repository root). Explicit entries are: `(0)`, `(1,m)`, `(2,J,m)`,
`(3,J13,J,m)`. Repeated components are deduplicated. J is derivative spin,
not total meson spin. Physical signed m is used; Redstar row = J-m+1.
All orders 0--3 produce 40 components, or 26 with negative m dropped.

FileDB retains `(t_slice, displacement, mom)` binary key layout; displacement
now holds an explicit coupled descriptor. Values have type 12, as in the
Redstar YLM prototype. FileDB output must be a new filename, to avoid collisions
with Cartesian keys and stale metadata. Superb storage records the same explicit
list along its d dimension. Metadata identifies `cg_coupled_covariant_derivatives_v1`, encoding
version 1, and the reduction flag. Existing Cartesian readers need an update
before consuming this new basis; changing only the producer is not enough.

Derivatives use the Redstar circular i factors, standard real CGs, and
`(1 x 1)->J` or `(1 x 1)->J13; (J13 x 1)->J`. Tensor positions 1 and 3 are
coupled first for three derivatives. Paths are stored in application order.
There is no extra factor -1/2 per derivative. Spin, subduction, and source/sink
phases are not included in these derivative-only matrices.

The existing superb prefix tree applies unit-link gauge derivatives to full
coordinate-space x vector-column tensors. It shares coordinate-space prefixes.
The callback accumulates the small Cartesian elemental matrices into coupled
components in double precision and writes each as soon as complete. This first
implementation still contracts the required Cartesian paths; it does not yet
couple coordinate-space tensors before contraction. Memory for accumulators is
bounded by the configured time/momentum block and requested components.

## Negative-m reconstruction

For equal left and right vector phasings and unit displacements:

    Phi[n,alpha,J,m](p)^dagger = eta * (-1)^m * Phi[n,alpha,J,-m](-p)
    eta = 1                       (n = 0, 1, 2)
    eta = (-1)^(3-J+J13)          (n = 3)

Thus the reader reconstructs negative m from positive -m at opposite momentum,
using a matrix adjoint, not elementwise conjugation. At p=0 this uses the same
momentum. With unequal vector phasings, the adjoint swaps those phasings as well.
The initial `drop_negative_m=true` implementation consequently requires equal
vector phasings and a momentum list containing every exact negative momentum.
It intentionally rejects asymmetric momentum lists even when periodic momentum
aliases could be used. The default is false. m=0 is retained for all momenta.

## Reproducible checks

From the Chroma repository root:

```sh
c++ -std=c++17 -Wall -Wextra -pedantic -Ilib tests/meson_ylm/dump_basis.cc \
  lib/meas/inline/hadron/meson_derivative_ylm.cc \
  lib/util/ferm/clebsch.cc -o /tmp/dump_meson_ylm
/tmp/dump_meson_ylm > /tmp/meson_ylm_basis.txt
python3 tests/meson_ylm/check_basis.py /tmp/meson_ylm_basis.txt /path/to/corr_graph.inspect.xml
```

Python requires numpy; omit the last argument to run without the external
Redstar fixture. Checks include descriptor validation, orthonormality of all
40 components, ordered adjoint identities, rectangular-column calculations on
noncommuting random SU(3) links, finite-momentum integration by parts, and
negative-m reconstruction at zero/nonzero momenta and nonzero equal vector
phasings. The optional XML check covers every vertex and source/sink pair,
including spin entries and reconstruction in the C++ derivative basis. It fits
subduction weights; it is not an independent subduction-phase calculation.
These tests do not replace a full Chroma gauge/colorvec run and storage readback.

## Current validation record

The full and reduced nonzero-momentum Chroma S3T runs now pass readback,
checksum, Cartesian reconstruction, and negative-m reconstruction checks.
See [validation](../../../hadron_elementals/notes/VALIDATION.md) and the [numerical archive](../../../hadron_elementals/numerics/README.md).
FileDB runtime and Redstar reader integration remain outstanding.

The parser still accepts singleton order-expansion requests such as `<elem>3</elem>`.
These are shorthand for all components at that derivative order, not explicit
coupled descriptors. The examples use explicit tuples to match the tested inputs.

The meson helper is separately compiled in `meson_derivative_ylm.cc` and uses
the same doubled-spin `Hadron::clebsch` implementation as the baryon helper.
The XML labels and meson circular/ordering conventions are unchanged.

## Storage precision (27 September)

S3T now defaults to 32-bit payloads; request `storage_precision=64` for references.
FileDB remains 64-bit. The existing runtime inputs explicitly retain 64; new
`.single.ini.xml` variants select 32 and distinct output names. See
[storage validation](../ylm_storage/README.md).
