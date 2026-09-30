## Regression tests for the propaneStackWarpInVar() correctness fixes.
##
## All three bugs only appear across a batch boundary or under a specific
## option combination, so each test compares against a configuration whose
## answer is known-good rather than pinning absolute pixel values.

.fr = function(nax = 900L) {
  suppressMessages({library(Rfits); library(Rwcs)})
  d = system.file('extdata/stack', package = 'ProPane')
  list(il = Rfits_make_list(dirlist = d, extlist = 2),
       kv = Rwcs_setkeyvalues(CRVAL1 = 36.8962, CRVAL2 = -5.1906, pixscale = 0.3,
                              NAXIS1 = nax, NAXIS2 = nax))
}

.stack = function(f, nb = 4, ...) suppressMessages(propaneStackWarpInVar(
  image_list = f$il, keyvalues_out = f$kv, inVar_list = 0.000468, cores = 2,
  multitype = 'fork', magzero_in = 30, magzero_out = 23.9, Nbatch = nb, ...))

test_that("weight_list is not indexed with the batch-local counter", {
  skip_on_cran()
  f = .fr()
  mat = lapply(5:8, function(i) matrix(3, dim(f$il[[i]])[1], dim(f$il[[i]])[2]))
  cases = list(
    'differing scalars'   = list(1, 2, 3, 1, 2, 3, 1, 2),
    'mixed scalar+matrix' = c(rep(list(2L), 4), mat)
  )
  for (nm in names(cases)) {
    many = .stack(f, nb = 8, weight_list = cases[[nm]])
    few  = .stack(f, nb = 4, weight_list = cases[[nm]])
    expect_identical(as.numeric(many$weight$imDat), as.numeric(few$weight$imDat),
                     info = nm)
  }
})

test_that("shapes unaffected by the bug stay unaffected", {
  skip_on_cran()
  f = .fr()
  mkmat = function(vals) lapply(seq_along(vals), function(i)
    matrix(vals[i], dim(f$il[[i]])[1], dim(f$il[[i]])[2]))
  for (wl in list(2L, mkmat(rep(2, 8)), mkmat(c(1, 2, 3, 1, 2, 3, 1, 2)))) {
    a = .stack(f, nb = 8, weight_list = wl)
    b = .stack(f, nb = 4, weight_list = wl)
    expect_identical(as.numeric(a$weight$imDat), as.numeric(b$weight$imDat))
  }
})

test_that("keep_extreme_pix + doclip resets hot before the restack", {
  skip_on_cran()
  f = .fr(nax = 2200L)
  noclip = .stack(f, nb = 8, keep_extreme_pix = TRUE)
  clip   = .stack(f, nb = 8, keep_extreme_pix = TRUE, doclip = TRUE)
  cp = as.numeric(clip$clip$imDat)
  hn = as.numeric(noclip$hot$imDat); hc = as.numeric(clip$hot$imDat)
  cn = as.numeric(noclip$cold$imDat); cc = as.numeric(clip$cold$imDat)
  at = !is.na(hn) & !is.na(hc) & cp > 0
  expect_gt(sum(at), 0L)                      # the clip path must have run
  # Clipping removes a frame, so hot can only fall and cold can only rise.
  expect_true(all(hc[at] <= hn[at] + 1e-9))
  expect_true(all(cc[at] >= cn[at] - 1e-9))
  # And the fix has to be visible: cold already behaved, hot must now too.
  expect_gt(sum(hc[at] < hn[at] - 1e-9), 0L)
})

test_that("multitype = 'cluster' accepts a cluster object", {
  skip_on_cran()
  suppressMessages(library(parallel))
  f = .fr()
  cl = makeCluster(2L)
  on.exit(stopCluster(cl), add = TRUE)
  s = suppressMessages(propaneStackWarpInVar(image_list = f$il,
        keyvalues_out = f$kv, inVar_list = 0.000468, cores = cl,
        multitype = 'cluster', magzero_in = 30, magzero_out = 23.9, Nbatch = 8))
  expect_gt(sum(as.numeric(s$weight$imDat), na.rm = TRUE), 0)
  ref = .stack(f, nb = 8)
  expect_identical(as.numeric(s$image$imDat), as.numeric(ref$image$imDat))
})

test_that("an integer warpgrid takes the coarse path", {
  skip_on_cran()
  f = .fr()
  src = suppressMessages(Rfits_read_image(file.path(
    system.file('extdata/stack', package = 'ProPane'), 'image_1.fits')))
  ex = suppressMessages(propaneWarp(src, keyvalues_out = f$kv, dim_out = c(900, 900),
       warpfield_return = TRUE, cores = 1))
  ap = suppressMessages(propaneWarp(src, keyvalues_out = f$kv, dim_out = c(900, 900),
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

## weight_image is derived from weight_list, so it must be trimmed by
## which_overlap in step with it. Left untrimmed it indexes frames by input
## position while weight_list indexes by trimmed position, so from the first
## dropped frame the two disagree.

test_that("weight_image is trimmed with weight_list across the overlap cut", {
  skip_on_cran()
  suppressMessages({library(Rfits); library(Rwcs)})
  d = system.file('extdata/stack', package = 'ProPane')
  # Read frames in rather than pointing at them: pointers carry per-read state
  # (XCUT*) that makes repeated runs in one session non-reproducible.
  read_all = function(shift4 = FALSE) {
    il = lapply(file.path(d, paste0('image_', 1:8, '.fits')),
                function(f) suppressMessages(Rfits_read_image(f)))
    if (shift4) {
      il[[4]]$keyvalues$CRVAL1 = il[[4]]$keyvalues$CRVAL1 + 5
      il[[4]]$keyvalues$CRVAL2 = il[[4]]$keyvalues$CRVAL2 + 5
    }
    il
  }
  kv = Rwcs_setkeyvalues(CRVAL1 = 36.8962, CRVAL2 = -5.1906, pixscale = 0.3,
                         NAXIS1 = 900, NAXIS2 = 900)
  run = function(iml, wll) suppressMessages(suppressWarnings(propaneStackWarpInVar(
    image_list = iml, keyvalues_out = kv, weight_list = wll, cores = 1,
    multitype = 'fork', magzero_in = 30, magzero_out = 23.9, Nbatch = 2)))

  # A matrix weight on ORIGINAL frame 5, which becomes trimmed index 4 once
  # frame 4 falls out of the overlap; every other frame gets a scalar. The
  # matrix must match frame 5's own dimensions (the frames are not uniform).
  d5 = dim(read_all()[[5]]$imDat)
  mk5 = matrix(1, d5[1], d5[2])
  wl = lapply(1:8, function(i) if (i == 5) mk5 else 2L)

  keep = c(1L, 2L, 3L, 5L, 6L, 7L, 8L)
  auto = run(read_all(TRUE), wl)
  expect_identical(as.integer(auto$which_overlap), keep)
  expect_equal(auto$Nim, 7L)
  # Hand-trimming is unambiguously aligned, so it is the known-good reference.
  hand = run(read_all(TRUE)[keep], wl[keep])
  expect_identical(as.numeric(auto$image$imDat), as.numeric(hand$image$imDat))
  expect_identical(as.numeric(auto$weight$imDat), as.numeric(hand$weight$imDat))
})

## The else branch indexed the whole mask_list rather than element i, so a
## plain-matrix mask hit `'list' object cannot be coerced to type 'double'`.

test_that("a plain-matrix mask_list is applied", {
  skip_on_cran()
  suppressMessages({library(Rfits); library(Rwcs)})
  d = system.file('extdata/stack', package = 'ProPane')
  il = lapply(file.path(d, paste0('image_', 1:8, '.fits')),
              function(f) suppressMessages(Rfits_read_image(f)))
  kv = Rwcs_setkeyvalues(CRVAL1 = 36.8962, CRVAL2 = -5.1906, pixscale = 0.3,
                         NAXIS1 = 900, NAXIS2 = 900)
  # Per-frame dimensions are not uniform, so each mask must match its own frame.
  mk = lapply(il, function(x) matrix(0L, dim(x$imDat)[1], dim(x$imDat)[2]))
  # An all-zero mask must leave the stack unchanged.
  base = suppressMessages(suppressWarnings(propaneStackWarpInVar(
    image_list = il, keyvalues_out = kv, cores = 1, multitype = 'fork',
    magzero_in = 30, magzero_out = 23.9, Nbatch = 8)))
  masked = suppressMessages(suppressWarnings(propaneStackWarpInVar(
    image_list = il, keyvalues_out = kv, mask_list = mk, cores = 1,
    multitype = 'fork', magzero_in = 30, magzero_out = 23.9, Nbatch = 8)))
  expect_identical(as.numeric(masked$image$imDat), as.numeric(base$image$imDat))
  # And a fully-masked frame must contribute nothing.
  mk_all = lapply(il, function(x) matrix(1L, dim(x$imDat)[1], dim(x$imDat)[2]))
  full = suppressMessages(suppressWarnings(propaneStackWarpInVar(
    image_list = il, keyvalues_out = kv, mask_list = mk_all, cores = 1,
    multitype = 'fork', magzero_in = 30, magzero_out = 23.9, Nbatch = 8)))
  expect_true(all(is.na(as.numeric(full$image$imDat))))
  expect_equal(sum(as.numeric(full$weight$imDat), na.rm = TRUE), 0)
})

## With multitype = 'cluster' the Nbatch default inherits the cluster object, so
## Nbatch was a list and every later comparison died before the projection.

test_that("multitype = 'cluster' works with the default Nbatch", {
  skip_on_cran()
  suppressMessages({library(Rfits); library(Rwcs); library(parallel)})
  f = .fr()
  cl = makeCluster(2L)
  on.exit(stopCluster(cl), add = TRUE)
  s = suppressMessages(propaneStackWarpInVar(image_list = f$il,
        keyvalues_out = f$kv, inVar_list = 0.000468, cores = cl,
        multitype = 'cluster', magzero_in = 30, magzero_out = 23.9))
  expect_equal(s$Nim, 8L)
  expect_gt(sum(as.numeric(s$weight$imDat), na.rm = TRUE), 0)
})
