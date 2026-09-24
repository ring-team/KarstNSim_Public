# Validation on 2026-09-24

Linux, GCC 13.3, CMake, .NET SDK 10.0.110. This is a tested local branch,
not an upstream release. Windows/macOS builds have not been verified.

## Runtime and compatibility evidence

- Release: all 11 CTest entries passed. After the last regional review fixes,
  all five SDK entries passed again.
- AddressSanitizer + UndefinedBehaviorSanitizer: all 11 entries passed with leak
  detection and no disabled bool checks. All five SDK entries also passed after
  the final regional changes. No sanitizer findings remained.
- Scientific regression: eight scenarios and 33 exports byte-identical to the
  optimization baseline. The in-memory API reproduced nine karst exports and
  eight connectivity matrices from the same corpus.
- Native concurrency: 14 checks under ThreadSanitizer, including independent
  and nested jobs. This machine required ASLR disabled for that TSan run.
- Installed C and C++ consumers passed against both static and shared builds.
  A .NET program called the shared C ABI and generated a region. The shared
  install was relocated to a path containing spaces; its CLI and consumers ran
  there. Vendored JSON symbols are hidden from the dynamic export table.
- Unmodified MOOCoW commit `eb3318a4712c2f360e8608c4167f4b84f062a294` accepted
  generated and constructed CSV fixtures with zero validator issues. Import,
  CSV re-export, SQLite save/reload and stored revision hashes matched. The
  existing seed map, schema and migration history stayed unchanged. Duplicate
  record IDs caused a complete transaction rollback. All DB writes were in a
  disposable test database; the actual MOOCoW checkout was not modified.

## Measured memory

Full-process peak RSS; decimal MB. Times are single observations on this
machine, not performance guarantees or controlled statistical comparisons.

| Workload | Peak RSS | Wall time | Result |
| --- | ---: | ---: | --- |
| Scientific full base, 549,999 support points | 722.57 MB | 20.04 s | Four exports match the prior full configuration |
| Coarse game preset, 24,278 support points | 44.97 MB | 0.99 s | Four exports match the same coarse preset |
| 64 sequential default regions, retaining only current/previous | 28.11 MB | 12.18 s | 768 caverns, 831 tunnels including stitches, 92.58 km of local passage |

The full scientific workload meets the 1 GB target and the named game workloads
meet the 100 MB target. The coarse preset changes sampling resolution relative
to the full base. Work/point/edge limits do not impose a universal allocator
cap. A whole world containing millions of miles was not generated here.

## Review corrections

- Native work consumption is reported on return/exception and charged to one
  regional work budget. Exact-budget and one-unit-under tests enforce it.
- Validation is cancellable during generation and bounded for public calls.
  Its overlap scan does not enumerate irrelevant authored/authored pairs.
- Generated facings survive MOOCoW's normalization exactly; authored values
  that would change are rejected. UTF-8 and storage text constraints are checked.
- Native geological scalar fields are initialized before being copied.
- Support-window indices include lattice origins; seam margins and port
  placement retries handle the tested boundary configurations.
- C ABI serialization avoids a second expensive geometry validation. JSON
  failures use the public error type. Completion manifests use atomic rename.
- Scientific iteration counts do not trigger an up-front result reservation.
  Shared libraries have explicit versions and relative runtime search paths;
  diagnostic compiler flags do not leak into installed CMake exports.

Detailed reproducible commands are in the SDK README, the adapter guide and
the companion handoff. The handoff records source/binary hashes, raw benchmark
metrics, logs and review reports. No rendered 3D wall mesh, application scheduler,
database merge policy, collision bake or navigation integration is claimed.
