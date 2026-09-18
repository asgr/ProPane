# ProPane 1.10.2

## Speed

`propaneWarp()` is substantially faster. Both changes are opt-in or bit-identical;
the default output is unchanged.

- `propaneWarpProPane()` now builds the warp field once and reuses it across all
  bands of a detection (new `warpfield_share = TRUE`, on by default). Previously
  each of the up-to-7 bands rebuilt an identical field. Output is bit-for-bit
  identical to 1.10.1. On a 5-band test detection this is a ~2.2x speedup, and
  ~2.4x over an 8-frame stack.

- New `warpgrid` argument to `propaneWarp()`, defaulting to `"exact"` which
  preserves 1.10.1 behaviour. `warpgrid = "auto"` (or a numeric lattice step)
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
