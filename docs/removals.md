# Reorg removal list

These removals apply to `rge/devel-datalib-reorg`.

| Component | Status | Decision |
| --- | --- | --- |
| `STOCH_HADRON` (`InlineStochHadronEnv`) | Removed, 9 November 2026 | User-directed removal of the legacy noise-dilution meson/baryon operator measurement. Its source, public header, registration and Autotools/CMake entries are deleted. No compatibility implementation, new storage adapter or format migration is provided for it. |

Inputs requesting `STOCH_HADRON` are no longer accepted by Chroma's measurement
factory. Historical ChangeLog entries remain as history. Other stochastic
measurements and the current distillation/SUPERB producers are separate code.
