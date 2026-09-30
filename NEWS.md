# ProPane 1.12.0

## Correctness

- `propaneStackWarpInVar()`: the stacking loops indexed `weight_image` and
  `weight_list` with the *batch-local* counter instead of the global frame index,
  so with `Nbatch < Nim` batch 2 reused batch 1's classification and read the
  wrong `weight_list` element. Which `weight_list` shapes are affected:

  | `weight_list`                 | broken at `Nbatch < Nim`? |
  | ----------------------------- | ------------------------- |
  | mixed scalar / matrix entries | **yes**                   |
  | differing scalar values       | **yes**                   |
  | one scalar for every frame    | no                        |
  | all matrices (same or differ) | no                        |

  All-matrix lists are safe because `pre_stack_weight_list` is rebuilt per batch
  and so is correctly indexed batch-locally, and the scalar/matrix test lands on
  a `TRUE` entry whichever index is used; a single uniform scalar is safe because
  every index carries the same value. Confirmed by measurement, not by reading.
  The image was unaffected in all cases and only the weight layer was wrong, so
  this silently produced a bad weight map. Four loops fixed, and the pixel-mask
  dumps on the re-project path now use the global index too (they were
  overwriting each other across batches).

- `propaneStackWarpInVar()`: `post_stack_hot[] -Inf` in the clip reset block
  computed `-Inf` and discarded it (there is no assignment), so with
  `doclip = TRUE` and `keep_extreme_pix = TRUE` the restack compared against the
  first pass's extremes, could never register a new hot pixel, and left clipped
  values in the returned `hot` layer. Now `post_stack_hot[] = -Inf`. The `hot`
  output changes for that combination: 54 of 4.8M pixels on the test frames.

- `propaneStackWarpInVar()` documented `multitype = 'cluster'` and accepted a
  cluster object for `cores`, but then ran `if(cores > Nbatch)` on it. A cluster
  is a list, so the comparison died with `'list' object cannot be coerced to type
  'integer'` for *every* call on that path, at any `Nbatch` and any cluster size
  -- the documented option had never worked. The guard now applies only when
  `cores` is a count; a cluster the caller built is left at the size they chose,
  since silently truncating the user's own object would be worse. Verified: the
  cluster path now runs and returns output bit-identical to `fork` for `image`,
  `weight`, `inVar` and `exp`. Note `propaneStackWarpMed()` and
  `propaneStackWarpQuanWt()` still `assertIntegerish(cores)` and so reject a
  cluster object outright; that looks deliberate and is left alone here.

  That fix was incomplete: the `Nbatch` default is `cores`, so a cluster caller
  who left `Nbatch` alone still died on the very next line (`if(Nbatch > Nim)` on
  a list). The guard only ever helped callers who passed an explicit numeric
  `Nbatch`. `Nbatch` now takes the cluster length when it inherits a non-numeric
  `cores`. The progress message also printed the cluster object rather than its
  size.

  Left alone deliberately: `cores` is shrunk to `Nbatch` *after*
  `registerDoParallel()` has already run, so that shrink affects only the message
  and never the worker count. Moving the registration below the sizing would make
  the shrink real, but it is a behaviour change with a performance trade-off of
  its own (with `return_all = TRUE`, `Nbatch` is raised to `Nim` afterwards, so
  the shrink would bite against a batch size that no longer applies). It is
  pre-existing and orthogonal to the crash this commit was fixing.

- `propaneStackWarpInVar()`: `weight_image` was derived from `weight_list` but
  was not trimmed by `which_overlap` alongside it. Once any frame failed the
  overlap check the two vectors disagreed on which frame was which, so from the
  first dropped frame a scalar weight could be read as a matrix or vice versa --
  either a silently wrong weight layer or a hard
  `pre_weight has non-conforming size` error. This is the same class of
  index-alignment bug as the batch-local counter fix above, on the other axis
  (input position vs trimmed position rather than global vs batch-local).

- `propaneStackWarpInVar()`: the `mask_list` fallback branch read `mask_list != 0`
  instead of `mask_list[[i]] != 0`, so a documented plain-matrix mask died with
  `'list' object cannot be coerced to type 'double'`. Present since 1.0.0, and
  now covered on both the single-pass and re-project paths.

- `.warpfield_coarse()` refused tolerances it could actually have met. The
  cost-floor check ran on the *SAFETY-adjusted* predicted step, which aims
  deliberately finer than the `s^2` law implies, so the prediction could fall
  below the floor while the floor itself still satisfied the request: on a
  900x900 field `warptol = 2.5e-8` predicted a sub-floor step from `s = 8` and
  was refused, yet step 4 measures 1.78e-08 and passes. The floor is now tested
  against the un-safetied prediction, and is clamped to rather than refused at,
  which restores the design intent that the prediction only chooses where to look
  and every candidate is still verified. A tolerance the floor genuinely cannot
  reach is still refused, and now on the first probe rather than after building
  a floor-sized lattice that would have failed anyway.

## Speed

`propaneStackWarpInVar()` is faster and uses less memory. All comparisons below
use the bundled 8-frame stack (1816x1836 in, 2200x2200 out) at `cores = 4`,
timed through an identical `pkgload` harness so both arms are compiled the same
way (columns are elapsed / user / sys seconds). The default `warpgrid = "exact"`
output is bit-for-bit unchanged: 53 layer comparisons against the previous build
were identical, and the only differences were the two deliberate correctness
fixes above. As a control, `cores = 1` is unchanged (10.77s -> 10.66s), which is
expected -- with one worker nothing crosses a pipe, so everything gained below is
inter-process communication.

- The image, inverse-variance, exposure and weight projections now happen inside
  a **single** parallel pass instead of four. The warp field for a frame is
  ~`nx * ny * 2 * 8` bytes (65 MB for a 2060x2070 field) and is identical for all
  four projections, but it used to be returned to the master (alongside every
  other frame's) and then sent back out to each worker as an argument. Holding it
  in a worker-local variable keeps it off the wire entirely: 5.67 / 19.9 / 13.1
  becomes 5.32 / 19.5 / 9.9, and peak RSS falls 2.62 GB -> 2.21 GB. System time is
  where the cost lived. At `cores = 8` the gain is larger (4.84 -> 4.28s, sys
  11.4s -> 6.3s, RSS 2.19 -> 1.64 GB) because more workers pay the round trip per
  frame. `Nbatch` also stops mattering for speed: the old code was 5-7% slower at
  `Nbatch = 4` than at 8, the new code is within 2%, so `Nbatch` is now
  essentially a memory knob rather than a speed knob too.

- `warpgrid` accepts an integer as its documentation always claimed. The branch
  was `if(warpgrid == 'approx')`, which is `FALSE` for every numeric value, so
  `warpgrid = 64` silently took the exact path and the fixed-step option did
  nothing.

- The default `warptol` is now `1e-3` rather than `1e-5`. This matters because at
  `1e-5` the coarse path could *never* be a win on a realistic field: bilinear
  error falls as the square of the lattice step, so meeting `1e-5` needs a step of
  ~4, and a step-4 field costs as much to build as the exact field while still
  being an approximation. It is still *accepted* at that tolerance (it resolves at
  the step-4 floor and honours its error claim), just no longer faster than
  `"exact"`, so the default is the thing that makes the coarse path worth using. Measured on a 2060x2070 field, exact build is 0.70s
  against 0.58s at step 4 and 2.35s at step 2. The old code discovered this by
  halving 64 -> 32 -> 16 -> 8 -> 4 and paying for all of it: `warpgrid = 'approx'`
  at the default tolerance was *slower* than `"exact"` (3.85s vs 3.65s end to
  end). `1e-3` resolves at step 29 in 0.033s against 0.70s for the exact field.
  Only reached when `warpgrid` is not `"exact"`, so the default output is
  untouched.

- `.warpfield_coarse()` predicts the lattice step it needs from one probe instead
  of halving blindly, and refuses rather than outspending the exact field. It
  also distinguishes *why* it gave up: 'the requested tolerance needs a lattice
  finer than step 4' is actionable, 'did not converge' is not. End to end with
  `warpgrid = 'approx'` end to end: 6.00s -> 3.95s at `cores = 4` (sys
  13.2 -> 10.8, RSS 2.61 -> 1.90 GB). Combining both changes the stack is 1.4x
  faster than the old default at `cores = 4` (5.67 -> 3.95s), and approx is now
  faster than exact rather than slower than it.

- `propaneStackWarpMed()` / `propaneStackWarpFunc()`: three `foreach` loops that
  ran *inside* each chunk task (the per-chunk overlap test, the frame cut, and
  the weight build) are now plain `vapply`/`lapply`. A nested `foreach` pays a
  fork/scheduler round trip that dwarfs the handful of `Rwcs_overlap` calls it
  wraps. Output is bit-for-bit identical across 12 configurations (compressed,
  uncompressed and partially-overlapping inputs at `chunk` 1000/500/220/110).
  Modest but free: uncompressed `chunk = 110` 9.44s -> 8.35s (11%),
  `chunk = 500` 1.45s -> 1.22s (16%), compressed `chunk = 500` 2.54s -> 2.29s.

## Output changes

- `propaneStackWarpInVar(return_all = TRUE)` no longer attaches the warp field to
  each projected pre-stack frame. The field is a by-product of projection, not a
  result, and at ~65 MB per frame it dominated the returned object: on the bundled
  frames `image_pre_stack` is 262 MB where it was 785 MB. Pixel values, headers and
  every other returned layer are unchanged, and nothing in the package or the
  vignettes reads `warpfield` off a pre-stack frame. If you do use it, take the
  field from `propaneWarp(..., warpfield_return = TRUE)` directly.

## Environments

- The warp and the WCS maths in `Rwcs`/`imager` are OpenMP-parallel, so each
  worker spawns a thread pool of its own. With `multitype = 'fork'` and
  `cores = 4` on a 16-core machine that is 4 workers x 16 threads doing very
  small jobs, and the thread churn dominates: a `cores = 1` stack burned 57.9s of
  CPU for 12.2s of wall time, against 12.2s CPU for 12.2s wall with
  `OMP_NUM_THREADS = 1`. Setting `OMP_NUM_THREADS` to `1` (or to `ncores / cores`)
  and letting `cores` do the parallelism is close to free and can be a large
  win; `cores = 8` with `OMP_NUM_THREADS = 1` beat `cores = 4` with the default
  on both wall time and total CPU. This is deliberately not set in the package:
  a user running `propaneWarp()` alone on one large frame wants those threads.

## Known limitations

- The dominant remaining cost of `propaneStackWarpMed()` scales with the *number
  of chunks*, i.e. `ceil(NAXIS1/chunk) * ceil(NAXIS2/chunk)`, not with the pixels
  that end up in the stack: on the same frames, `chunk = 2200` takes 0.74s and
  `chunk = 137` takes 6.30s for identical output. With compressed inputs every
  windowed read decompresses the rows it touches regardless of window width (a
  100x100 window costs 0.018s against 0.031s for the whole frame), so small
  chunks are especially expensive. Grouping horizontally adjacent chunks so a
  task reads one wide window and slices it in memory was measured and is
  numerically safe (a slice of a wide read is bit-identical to the narrow read,
  NA edge padding included) but is *not* implemented, because for uncompressed
  pre-dumped frames reads are only 2-11% of the loop and the multiplication is
  read-independent. Users who need small `chunk` values should watch this space.

## Infrastructure

- Added `bench/stack-one.R`, `bench/ab-sweep.sh` and `bench/make-variants.sh`,
  which build four variants of the package (pristine HEAD, warp-only, stack-only,
  both) into `/tmp/ppvar` so the two changes can be attributed separately under an
  identical `pkgload` harness. Comparing against the *installed* build is
  misleading: it is byte-compiled and `pkgload` is not, which alone accounts for
  ~25% and initially inverted the sign of the fusion result.
- Added `bench/identity-check.R`, which A/Bs every returned layer of eight
  stacking configurations against the installed build and labels the differences
  that are deliberate correctness fixes.
- Added `bench/approx-accuracy.R`, which reports flux conservation, per-pixel
  delta in units of the sky rms, NA-pattern and weight changes for each
  `warptol`.
- Added `bench/stackbench.R` for the whole-stack time/CPU/RSS sweep.

# ProPane 1.11.0

## Correctness

- When the requested output frame does not overlap the input, `propaneWarp()`
  now returns the full requested frame filled with `blank` and says so, instead
  of warping. The tight-crop range for a non-overlapping frame is inverted
  (`x [9998, 1816]` over an 1816-wide input), and that pair was being passed to
  a *box* crop, which read it as a width and silently expanded the working grid:
  a shift of 10000 px in `CRPIX1` built an 8183x1768 = 14.5M pixel warp field
  (9.7s, much of it inside wcslib non-convergence) against a 3.2M pixel
  reference, then returned a normal-looking block of NaN. Note this changes the
  output of such calls: it is now `blank` rather than NaN, which also means
  `blank` is finally honoured (it was previously overwritten by the warp), and
  the returned header keeps the requested `CRPIX` and a full-frame crop rather
  than inheriting the expansion. A warp with any overlap is unaffected.

## Speed

`propaneWarp()` is substantially faster. Both changes are opt-in or bit-identical;
the default output is unchanged.

- `propaneWarpProPane()` now builds the warp field once and reuses it across all
  bands of a detection (new `warpfield_share = TRUE`, on by default). Previously
  each of the up-to-7 bands rebuilt an identical field. Output is bit-for-bit
  identical to 1.10.1. On a 5-band test detection this is a ~2.2x speedup, and
  ~2.4x over an 8-frame stack.

- New `warpgrid` argument to `propaneWarp()`, defaulting to `"exact"` which
  preserves 1.10.1 behaviour. `warpgrid = "approx"` (or a numeric lattice step)
  evaluates the celestial transform on a coarse lattice and reconstructs the
  field as a least-squares affine plus a bilinearly interpolated residual,
  rather than transforming every output pixel. The per-pixel transform is about
  three quarters of a call, so this is the large win: ~3x on a single band, and
  ~3.5x on an 8-frame 5-band stack (83.7s -> 24.3s).

  The reconstruction is verified against the true transform at nine interior
  points per lattice cell, and against a separate dense finiteness sweep; if the
  requested `warptol` cannot be met the lattice is refined and, failing that,
  the exact field is used with a message. A warp that refuses to be approximated
  costs time but is never less accurate than requested. On typical ASKAP-style
  fields the maximum warp-field error is ~5e-06 pixels and integrated flux agrees
  to ~1e-07 relative; individual pixels whose interpolation weights are truncated
  at the edge of the input can differ by more, so this is opt-in.

- A tight crop whose bounds invert (an output frame with little or no overlap
  with the input) previously passed an expanded, silently-growing working grid to
  the warp -- measured at 4.5x the pixel count and a normal-looking but all-empty
  result. It now warns. The returned data is unchanged.

## Known issue

`direction = "forward"` is not reproducible run-to-run in 1.10.1 and still is:
imager/CImg accumulates into the output without atomics, so results vary between
runs and depend on `OMP_NUM_THREADS`. Backward warping is unaffected. The golden
test for the forward case checks statistical agreement rather than bit-identity
to reflect this. Fixing it requires a change upstream in CImg.

## Infrastructure

- Added a `testthat` suite (`tests/`) for `propaneWarp()`, which previously had
  none. Output is pinned with raw-byte digests over 12 geometry/interpolation/
  boundary combinations, so future changes to this numerically-lossy workhorse
  are caught exactly rather than approximately.
- Added `bench/` with the warp benchmark, the golden-fixture generator, and a
  synthetic WCS stress gate for the coarse path.
- `digest` added to Suggests (tests only).
