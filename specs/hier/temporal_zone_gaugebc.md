# Temporal Zone GaugeBC Spec

## Goal

Add a new gauge boundary-condition class that uses the existing `GaugeBC::zero(P&)`
hook to zero gauge-like force/momentum fields on a configurable set of time
coordinate intervals.

Chosen public name:

- C++ class: `TemporalZoneGaugeBC`
- XML/factory name: `TEMPORAL_ZONE_GAUGEBC`
- Files: `lib/actions/gauge/gaugebcs/temporal_zone_gaugebc.h` and
  `lib/actions/gauge/gaugebcs/temporal_zone_gaugebc.cc`

This is a force-suppression boundary condition. It does not modify gauge links.

## Current Repo Facts

- `GaugeBC<P,Q>` is defined in `lib/gaugebc.h` and requires:
  - `modify(Q& u) const`
  - `zero(P& ds_u) const`
  - `nontrivialP() const`
- Gauge boundary conditions are constructed through
  `lib/actions/gauge/gaugebcs/gaugebc_factory.h`.
- The gauge BC registration aggregator is
  `lib/actions/gauge/gaugebcs/gaugebc_aggregate.cc`.
- Gauge actions already call `getGaugeBC().zero(...)` after force construction.
  Examples:
  - `lib/actions/gauge/gaugeacts/plaq_gaugeact.cc`
  - `lib/actions/gauge/gaugeacts/rect_gaugeact.cc`
  - `lib/actions/gauge/gaugeacts/aniso_sym_spatial_gaugeact.cc`
  - `lib/actions/gauge/gaugeacts/aniso_sym_temporal_gaugeact.cc`
- Existing mask-based zeroing is implemented by the Schroedinger classes:
  - `SchrGaugeBC::zero(...)` in `schroedinger_gaugebc.cc`
  - `SchrSFZeroGaugeBC::zero(...)` in `schr_sf_zero_gaugebc.cc`
  Both use a `multi1d<LatticeBoolean>` mask and `copymask`.

## XML Interface

Use this under any gauge action `<GaugeBC>` block:

```xml
<GaugeBC>
  <Name>TEMPORAL_ZONE_GAUGEBC</Name>
  <t_dir>3</t_dir>
  <zero_intervals>
    <elem>
      <t_start>0</t_start>
      <t_end>2</t_end>
    </elem>
    <elem>
      <t_start>7</t_start>
      <t_end>7</t_end>
    </elem>
  </zero_intervals>
</GaugeBC>
```

Semantics:

- `t_dir` is optional; default is `Nd - 1`.
- `zero_intervals` is required and contains one or more inclusive intervals.
- An interval selects all sites whose `Layout::latticeCoordinate(t_dir)` satisfies
  `t_start <= t <= t_end`.
- All link directions `mu = 0 ... Nd-1` are zeroed on selected sites.
- Intervals must satisfy `0 <= t_start <= t_end < Layout::lattSize()[t_dir]`.
- Wrapping intervals are rejected in the first implementation. Users should spell
  a wrapping interval as two ordinary intervals.

## Parameters

Add a small parameter type in the new header:

```c++
struct TemporalZoneInterval {
  int t_start;
  int t_end;
};

struct TemporalZoneGaugeBCParams {
  TemporalZoneGaugeBCParams();
  TemporalZoneGaugeBCParams(XMLReader& xml, const std::string& path);

  int t_dir;
  multi1d<TemporalZoneInterval> zero_intervals;
};
```

Implement `read/write` overloads for both structs so the `multi1d<...>` XML
reader can consume `<zero_intervals><elem>...</elem></zero_intervals>`.

## Class Behavior

The class should derive from:

```c++
GaugeBC<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix>>
```

It should own:

```c++
TemporalZoneGaugeBCParams param;
multi1d<LatticeBoolean> mask;
```

Constructor:

1. Store `param`.
2. Validate `t_dir` and each interval.
3. Build `mask.resize(Nd)`.
4. Build one site mask:
   ```c++
   LatticeInteger t = Layout::latticeCoordinate(param.t_dir);
   LatticeBoolean site_mask = false;
   for each interval:
     site_mask |= (t >= interval.t_start) && (t <= interval.t_end);
   ```
5. Assign `mask[mu] = site_mask` for every link direction `mu`.

Methods:

```c++
void modify(multi1d<LatticeColorMatrix>& u) const
{
  // No-op. This BC only suppresses force/update fields.
}

void zero(multi1d<LatticeColorMatrix>& ds_u) const
{
  LatticeColorMatrix z = QDP::zero;
  for (int mu = 0; mu < ds_u.size(); ++mu) {
    copymask(ds_u[mu], mask[mu], z);
  }
}

bool nontrivialP() const
{
  return param.zero_intervals.size() > 0;
}
```

## Registration

Add a namespace matching local style:

```c++
namespace TemporalZoneGaugeBCEnv {
  extern const std::string name;
  bool registerAll();
}
```

In `temporal_zone_gaugebc.cc`:

- Define `name = "TEMPORAL_ZONE_GAUGEBC"`.
- Register the create callback with `TheGaugeBCFactory`.
- Match the Schroedinger classes for the first implementation by registering the
  double-precision gauge field factory only. Add `TheGaugeBCFFactory` and
  `TheGaugeBCDFactory` later only if single/double-specific gauge action paths need it.

Update:

- `lib/actions/gauge/gaugebcs/gaugebc_aggregate.cc`
  - include the new header
  - call `TemporalZoneGaugeBCEnv::registerAll()`
- `lib/actions/gauge/gaugebcs/gaugebcs.h`
  - include the new header

## Build System Updates

Add the new header and source to both build lists:

- `lib/Makefile.am`
  - header list near the other `actions/gauge/gaugebcs/*.h`
  - source list near the other `actions/gauge/gaugebcs/*.cc`
- `lib/CMakeLists.txt`
  - header list near lines containing `schr_sf_zero_gaugebc.h`
  - source list near lines containing `schr_sf_zero_gaugebc.cc`

## Tests

Add at least one test input and one focused code test.

Suggested XML smoke input:

- Copy a small gauge-only leapfrog input such as
  `tests/t_leapfrog/t_leapfrog.plaq.ini.xml`.
- Replace `PERIODIC_GAUGEBC` with `TEMPORAL_ZONE_GAUGEBC`.
- Use a small lattice with `nrow` time extent at least 4 and one interval such
  as `[0,1]`.

Suggested focused test:

1. Instantiate `TemporalZoneGaugeBC` with `t_dir = Nd - 1` and intervals `[0,1]`.
2. Fill a `multi1d<LatticeColorMatrix> ds_u` with nonzero values.
3. Call `zero(ds_u)`.
4. Verify all directions are zero on `t = 0,1` and unchanged elsewhere.

If no standalone assertion utility is convenient, add diagnostic XML/log checks
using `sum(localNorm2(ds_u[mu]))` before and after masking.

## Validation Checklist

1. Build `chroma-hier`.
2. Verify the new factory name is accepted in a gauge action XML.
3. Verify force norms vanish on selected time intervals.
4. Verify force norms outside selected intervals are unchanged relative to the
   same input with `PERIODIC_GAUGEBC`.
5. Verify overlapping intervals behave like their union.
6. Verify invalid intervals abort with a clear message.

## Scope Notes

- This spec intentionally leaves `modify(Q&)` as a no-op, so it does not impose
  Dirichlet values on gauge links.
- Gauge-force paths get this behavior automatically because gauge actions call
  `getGaugeBC().zero(...)`.
- If the intended physics requires frozen link updates, not just zero force,
  check the trajectory/integrator path too. In the current `LatColMatHMCTrj` and
  `LatColMatSMDTrj`, refreshed momenta are projected with `taproj(...)` but are
  not passed through a gauge BC object. That is separate from this BC class and
  would need an additional trajectory or integrator hook.
