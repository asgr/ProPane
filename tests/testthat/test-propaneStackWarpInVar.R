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
