# AGENTS.md — working on ProPane

Guidance for AI agents and new contributors. Read this before editing.

## What this package is

ProPane warps and stacks astronomical FITS images. It is the successor to the
warp/stack functions that used to live in `Rwcs` and `ProFound`; `README.md` has
the old→new mapping table.

- **`propaneWarp()`** reprojects one WCS onto another, via `wcslib` (through
  `Rwcs`) for the celestial mapping and `CImg` (through `imager`) for pixel
  resampling. Flux-conserving, with photometric zero-point rescaling and
  ZIMAGE/tile-compressed FITS handling.
- **`propaneStackWarpInVar()`** is the flagship: a parallel
  inverse-variance-weighted stack of many dithered frames onto a common WCS,
  returning image / inVar / weight / exposure maps plus hot–cold pixel maps.
  Siblings: `propaneStackWarpMed()`, `propaneStackWarpQuanWt()`,
  `propaneStackWarpFunc()`, `propaneStackFlatInVar()`, `propaneWarpDump()`.
- **`propaneTweak()`** solves small astrometric offsets + rotations;
  `propaneWCSmod()` / `propaneTran()` apply the result.
- Utilities: `propaneFrameFinder()`, `propaneGenWCS()`, `propaneCovTest()`,
  `propaneBadPix()`/`propanePatchPix()`, `propaneDilate()`, `propaneLocalFunc()`,
  `propaneInterp2D()`/`propaneBin2D()`.

Heavy dependencies (`Rfits`, `Rwcs`, `imager`, `ProUtils`, `ProFound`, `Rfast2`)
are installed from GitHub (`asgr/*`), not CRAN — see `Remotes:` in `DESCRIPTION`.

## Layout

| Path | Notes |
|---|---|
| `R/` | one file per topic area (~17 files); `R/RcppExports.R` is generated, never edit |
| `src/` | 5 C++ files (`warpfield.cpp`, `stack.cpp`, `propaneInterp2D.cpp`, `dilate.cpp`, `mat_diff_sum.cpp`); no `Makevars`, no OpenMP |
| `man/` | **hand-written** `.Rd`, organised by *topic family*, not per function |
| `NAMESPACE` | **hand-maintained**; organised in three commented groups |
| `tests/testthat/` | golden-fixture regression tests; see `Testing` below |
| `inst/extdata/` | the 8 stacked FITS frames + `swarp_med.fits` — the only test/vignette data, do not modify |
| `bench/` | developer-only scripts, excluded from the build by `.Rbuildignore` |
| `vignettes/` | `ProPane_Stack.Rmd`, `ProPane_Tweak.Rmd` |

## This package does not use roxygen2

There are zero `#'` lines. Documentation is written by hand and `NAMESPACE` is
edited by hand. Do **not** run `devtools::document()`.

If you add or change an exported function's arguments, update all three by hand:

1. the definition in `R/`;
2. the `\usage{}` (and `\arguments{}`) block of the relevant family `.Rd` in
   `man/` — add a `\alias{}` there if the function is new to that family;
3. `NAMESPACE` if the function is newly exported, in the correct group comment.

Two known inconsistencies already in the tree: `propaneSegimShare` has a `.Rd`
but is not exported, and `man/propaneBadPix.Rd` / `man/propaneStackFlat.Rd`
alias functions (`propanePix`, `propaneStackFlat`) that no longer exist. Leave
them alone unless you are deliberately cleaning up.

## Code conventions

- **Use `=`, not `<-`.** The only exceptions in `R/` are the two `<<-` in the
  `foreach` loop in `R/propaneWarp.R` and generated `RcppExports.R`.
- **Exported functions are prefixed `propane`** (lowercase `p`, camelCase after).
  **Internal helpers are dot-prefixed** (`.warpfield_coarse`, `.cost_fn_image`,
  `.makeBrush`) and are never exported or documented.
- **Return `invisible()`** from exported functions that are called mostly for
  their effect or whose result is normally assigned.
- A `"ProPane"` object is a **plain list with a `class` attribute** — no S3
  methods, no `print`/`plot`/`summary`. Members include `image`, `weight`,
  `inVar`, `exp`, `cold`, `hot`, `clip`, `which_overlap`, `time`, `Nim`,
  `dump_dir`. Don't introduce a method-based class without asking.
- **Comment densely and explain *why*, with measurements where you have them.**
  This is the strongest convention in the codebase: the long block comments at
  the top of `R/propaneWarp.R`, `R/propaneStackWarpMed.R` and the test files are
  design notes, and they are what make the code reviewable. Match that style.
- British spelling in prose and messages, single-quoted strings.
- Validation is inconsistent: `checkmate::assert*` is used heavily in
  `propaneBadPix.R`, lightly elsewhere, and not at all in `propaneWarp()`. Follow
  the local file rather than trying to standardise the package.
- **Parallel work** goes through `foreach` + `doParallel`, registered per call
  with `registerDoParallel()` and always torn down with `stopImplicitCluster()`.
  Declare loop variables (`i = j = NULL`) so `R CMD check` doesn't warn.
- **Soft dependencies** (`imager`, `RANN`, `cmaes`, `ProUtils`, `Rfast2`) are
  called with `::` behind `if(!requireNamespace('x', quietly=TRUE)) stop(...)`.
  Never add them to `Imports` or to the blanket `import()` block.
- Pixel-centre convention is threaded as `pixcen = 'FITS'` (1-based) vs `'R'`
  (0.5 offset) in `propaneInterp2D()` / `propaneBin2D()`. Keep it explicit.

### C++ changes

Edit the `.cpp` file, then call `Rcpp::compileAttributes()` to regenerate
`R/RcppExports.R` and `src/RcppExports.cpp` — commit both. R wrappers are the
dot-prefixed functions in `R/RcppExports.R`; package-level wrappers with
validation live in the corresponding `R/propane*.R`.

## Testing

```r
devtools::test()                       # or testthat::test_local()
testthat::test_file("tests/testthat/test-propaneWarp.R")
```

`tests/testthat.R` is the standard three-liner, so `R CMD check` runs the same
suite. The CI matrix is macOS release + Ubuntu devel/release/oldrel-1; Windows is
deliberately not tested.

Tests are **golden-record regression tests**, not unit tests:

- `test-propaneWarp.R` compares outputs against
  `tests/testthat/fixtures/warp-golden.rds`: exact shape, exact NA/NaN pattern,
  sampled pixels to `GOLDEN_TOL = 1e-9` relative, and selected header keys.
- `test-warpgrid-stress.R` synthesises SIP/TPV/CAR/SIN headers in memory and
  asserts the coarse warp field either meets its tolerance or returns `NULL`.
- Digests use `tools::md5sum` on raw doubles, deliberately avoiding a `Suggests`
  package, because CI installs hard dependencies only.
- `direction = "forward"` is recorded as `nondeterministic` and gated loosely
  (flux + ≥50% of sampled pixels). CImg's forward scatter has no atomics, and
  whether the race is reachable depends on the CImg build: stock CRAN `imager`
  parallelises it, the `asgr/imager` fork serialises it. **Never decide this
  case's status by re-running it locally** — on a serial build every re-run
  agrees and the case will pass locally and fail on Ubuntu. That is precisely
  what happened in 1.12.1.

Rules when touching warp behaviour:

1. If your change is **not** meant to alter `propaneWarp()` output, the existing
   golden record must still pass unchanged. That is the point of the tolerance.
2. If your change **is** meant to alter output, regenerate the record *before*
   editing, from an installed build of the old code:
   `R CMD INSTALL . && Rscript bench/make-fixtures.R`, then commit the new
   fixture alongside the code change and say plainly in `NEWS.md` that pixels
   moved.
3. Keep the suite fast (currently ~35 s locally; CI is 4–5× slower). Prefer
   shrinking `dim_out` over dropping a case.
4. **Do not add `foreach`/`%dopar%` or `Rfits_make_list()` to tests.** A torn-down
   cluster backend leaves dead sockets that kill every later Rfits call in the
   session. This is why the end-to-end stacking tests were removed (see the
   header comment in `test-propaneStackWarpInVar.R`; the old tests are in git
   history at `5a345ee..e64e40c`).
5. **`skip_on_cran()` does not skip in this repo's CI** — r-lib/actions sets
   `NOT_CRAN=true`. Use an explicit guarded `skip()` if you need an out.

## Benchmarking

`bench/` is excluded from builds and only used by maintainers:

```sh
Rscript bench/warpbench.R        # warp timings; BENCH_REPS env var, default 7
Rscript bench/stress-warpgrid.R
Rscript bench/make-fixtures.R    # needs an *installed* build, not load_all()
```

Two scripts use `pkgload::load_all()`; `make-fixtures.R` uses
`library(ProPane)`. Report min and median timings, never a single run.

## Before you ship

```r
devtools::check(args = c("--no-manual", "--as-cran", "--no-build-vignettes"))
```

Checklist:

- [ ] `devtools::test()` passes and runtime hasn't ballooned.
- [ ] `NEWS.md` has a section for the new version, with entries under
      `## Correctness` / `## Speed` / `## Output changes` / `## Infrastructure`
      as applicable. Entries are expected to be specific and measurement-backed.
- [ ] `DESCRIPTION` `Version` and `Date` bumped together (the project bumps the
      patch version per change and keeps a separate `Version bump` commit).
- [ ] Manual `man/*.Rd` and `NAMESPACE` updated in step with any signature change.
- [ ] No new `Suggests` dependency for something the tests need at check time;
      CI installs hard dependencies plus a fixed `extra-packages` list.
- [ ] Vignettes still build (they need `imager` from GitHub).

## Commit style

Free-form, sentence-case imperatives: `Modified coarse field limit (4 -> 8)`,
`Some review fixes`, `Version bump`, `Fater warping using affine lattice
distortions.` Do **not** use conventional-commit prefixes (`feat:`, `fix:`).
Version bumps are their own commit. Substantive detail belongs in `NEWS.md`, not
in the commit body.

## Known rough edges

Don't be surprised by these; fix only if that is the task.

- `NEWS.md` lags `DESCRIPTION` (top heading is 1.12.0).
- `propaneCovTest()` has an unresolved `#something here!!!!!!` at `R/propaneCovTest.R:5`.
- `R/propaneCropOverlap.R` mixes `|` and `&` without parentheses in one condition.
- Dead code: ~70 lines of commented-out `propaneNifty` in `R/propaneTweak.R`, and
  commented-out `ProUtils`/`Rfast2` variants in the stacking files.
- Several `stop()` messages name the wrong missing package (copy-paste).
- `NEWS.md` mentions `bench/` scripts that no longer exist.