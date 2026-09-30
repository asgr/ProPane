## Tests for the coarse warp-field builder used by propaneStackWarpInVar().
##
## The end-to-end propaneStackWarpInVar() regression tests were removed: every
## one of them stacked the full eight-frame example set at 900^2 or 2200^2,
## several times per test, which cost ~140s locally and *hours* on the GitHub
## Actions Ubuntu runners. The Ubuntu jobs began hanging the moment those tests
## landed (commit 5a345ee) and have never completed since, while the macOS jobs
## running the identical suite finished in ~15min. All of them went through
## multitype = 'fork'/'cluster', so the fork-based doParallel backend is the
## prime suspect rather than raw cost.
##
## Note skip_on_cran() does NOT skip these under r-lib/actions, which sets
## NOT_CRAN=true -- they were running in full on every push.
##
## What survives here drives .warpfield_coarse() directly: no stacking, no
## fork, and about a second of runtime for the whole file. The removed tests
## are in git history at 5a345ee..e64e40c if they are wanted back behind a
## non-CI runner.

# Only the output keyvalues are needed now that the stacking tests are gone.
# Deliberately no Rfits_make_list(): it runs a foreach %dopar% loop, and if a
# doParallel cluster backend has been registered and then torn down by an
# earlier test, the sockets are dead and every later Rfits call in the session
# dies with "invalid connection". Keeping this file free of foreach means it
# cannot inherit that, nor hand it on.
.kv_out = function(nax = 900L) {
  suppressMessages(library(Rwcs))
  Rwcs_setkeyvalues(CRVAL1 = 36.8962, CRVAL2 = -5.1906, pixscale = 0.3,
                    NAXIS1 = nax, NAXIS2 = nax)
}

test_that("an integer warpgrid takes the coarse path", {
  skip_on_cran()
  suppressMessages(library(Rfits))
  kv = .kv_out()
  src = suppressMessages(Rfits_read_image(file.path(
    system.file('extdata/stack', package = 'ProPane'), 'image_1.fits')))
  ex = suppressMessages(propaneWarp(src, keyvalues_out = kv, dim_out = c(900, 900),
       warpfield_return = TRUE, cores = 1))
  ap = suppressMessages(propaneWarp(src, keyvalues_out = kv, dim_out = c(900, 900),
       warpfield_return = TRUE, cores = 1, warpgrid = 64L))
  # Before the fix `warpgrid = 64` fell through to the exact branch and was
  # bit-identical to the default; the coarse field must differ but agree closely.
  a = cbind(as.vector(ex$warpfield[, , 1]), as.vector(ex$warpfield[, , 2]))
  b = cbind(as.vector(ap$warpfield[, , 1]), as.vector(ap$warpfield[, , 2]))
  expect_identical(dim(a), dim(b))
  ok = is.finite(a) & is.finite(b)
  expect_false(identical(a[ok], b[ok]))
  expect_lt(max(abs(a[ok] - b[ok])), 1e-2)
})

test_that("a tolerance needing a sub-step-4 lattice refuses rather than grinding", {
  skip_on_cran()
  suppressMessages({library(Rfits); library(Rwcs)})
  ns = asNamespace('ProPane')
  coarse = get('.warpfield_coarse', envir = ns)
  out2in = get('.warpfunc_out2in', envir = ns)
  src = suppressMessages(Rfits_read_image(file.path(
    system.file('extdata/stack', package = 'ProPane'), 'image_1.fits')))
  ref = suppressMessages(Rfits_read_image(file.path(
    system.file('extdata/stack', package = 'ProPane'), 'image_2.fits')))
  wf = function(x, y, cores = 1, ...)
    out2in(x, y, header_in = Rfits_keyvalues_to_raw(src$keyvalues),
           header_out = Rfits_header_to_raw(Rfits_keyvalues_to_header(ref$keyvalues)),
           cores = cores)
  # 1e-9 px is unreachable at a step worth building, so refuse immediately.
  expect_null(coarse(wf, c(900, 900), tol = 1e-9, step0 = 64L, cores = 1))
  expect_null(coarse(wf, c(900, 900), tol = 1e-3, step0 = 2L, cores = 1))
  # A reachable tolerance still succeeds and honours its own error claim.
  r = coarse(wf, c(900, 900), tol = 1e-3, step0 = 64L, cores = 1)
  expect_false(is.null(r))
  expect_lt(r$maxerr, 1e-3)
  expect_gte(r$step, 4L)
})

## The cost floor clamps the predicted step rather than refusing on the
## prediction alone: predicting a sub-floor step only says the s^2 law
## extrapolates past it, not that step 4 would miss.

test_that("a tolerance met at the cost-floor step is not refused", {
  skip_on_cran()
  suppressMessages({library(Rfits); library(Rwcs)})
  ns = asNamespace('ProPane')
  coarse = get('.warpfield_coarse', envir = ns)
  out2in = get('.warpfunc_out2in', envir = ns)
  src = suppressMessages(Rfits_read_image(file.path(
    system.file('extdata/stack', package = 'ProPane'), 'image_1.fits')))
  ref = suppressMessages(Rfits_read_image(file.path(
    system.file('extdata/stack', package = 'ProPane'), 'image_2.fits')))
  wf = function(x, y, cores = 1, ...)
    out2in(x, y, header_in = Rfits_keyvalues_to_raw(src$keyvalues),
           header_out = Rfits_header_to_raw(Rfits_keyvalues_to_header(ref$keyvalues)),
           cores = cores)
  # At step 8 this tol predicts a lattice finer than the floor, yet step 4
  # measures 1.78e-08 and satisfies it. Refusing here would discard a field
  # that is both accurate enough and cheaper than exact.
  r = coarse(wf, c(900, 900), tol = 2.5e-8, step0 = 8L, cores = 1)
  expect_false(is.null(r))
  expect_equal(r$step, 4L)
  expect_lt(r$maxerr, 2.5e-8)
  # The floor is still respected: nothing finer than step 4 is ever built.
  expect_gte(r$step, 4L)
  # And a tolerance step 4 genuinely cannot meet is refused on the prediction,
  # without paying for a floor-lattice build that would fail anyway.
  expect_null(coarse(wf, c(900, 900), tol = 1e-8, step0 = 64L, cores = 1))
  expect_match(get('.warpfield_refusal', envir = ns)(), 'finer than the cost-floor step 4')
})
