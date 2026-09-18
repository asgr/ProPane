# Golden-fixture regression tests for propaneWarp().
#
# fixtures/warp-golden.rds was produced by released ProPane 1.10.1 (regenerate
# with bench/make-fixtures.R, run against the *unmodified* code). Steps 1-3 of
# the warp performance work must not change the output; only the opt-in
# `warpgrid` path is permitted to differ.
#
# "Must not change" is enforced to a tight tolerance rather than bit-for-bit,
# because exact digests do not survive a change of build environment. The
# residual between 1.10.1 and the current code built in a different Rfits/Rwcs
# environment is a *uniform* factor of ~8.4e-14, traceable to pixscale(): its CD
# matrix round-trips through 15-significant-digit FITS cards, and it is sampled
# at the field centre, where a 0.001 px shift already moves the scale by 2e-9.
# Verified neutral: the 1.10.1 source and HEAD, installed side by side in one
# environment, give byte-identical digests for every golden case.
#
# So a digest comparison would be asserting "this machine is unchanged" rather
# than "this code is unchanged". The gates below are instead: exact shape, exact
# missing-value pattern (NA count, NA positions, NA vs NaN), every sampled pixel
# agreeing to GOLDEN_TOL relative, and bit-stability across repeated runs.

fixture_path <- testthat::test_path("fixtures", "warp-golden.rds")

# Relative agreement required of the golden sample. Observed environment drift
# is ~1e-13; a genuine warp or crop bug moves pixels by O(1e-3) or breaks the NA
# pattern outright, so this sits four orders above the noise and six below the
# smallest regression worth catching.
GOLDEN_TOL <- 1e-9

.cache <- new.env(parent = emptyenv())

.golden <- function() {
  if (is.null(.cache$golden)) .cache$golden <- readRDS(fixture_path)
  .cache$golden
}

.frames <- function() {
  if (is.null(.cache$frames)) {
    suppressMessages(library(Rfits))
    suppressMessages(library(imager))
    d <- system.file("extdata/stack", package = "ProPane")
    g <- .golden()
    .cache$frames <- list(
      src = Rfits_read_image(file.path(d, g$src_frame)),
      ref = Rfits_read_image(file.path(d, g$ref_frame)),
      dim = g$dim_ref
    )
  }
  .cache$frames
}

# Same recipe as bench/make-fixtures.R: hash raw double bytes via tools::md5sum
# so the digest is independent of R's serialization format, sensitive to NA vs
# NaN, and available without any Suggests package (CI installs hard deps only).
.digest_matrix <- function(m) {
  v <- as.vector(m, mode = "double")
  f <- tempfile()
  con <- file(f, "wb")
  writeBin(v, con, size = 8)
  close(con)
  on.exit(unlink(f))
  unname(tools::md5sum(f))
}

# Compare a result against a golden record. Everything that must be exact is
# checked exactly; only the pixel values are allowed the tolerance. The +1 in
# the denominator keeps near-zero pixels (blank sky, and the magzero case where
# old values straddle 0) from inflating a relative error.
.expect_close_golden <- function(mat, g, label) {
  info <- paste(label, collapse = " ")

  expect_identical(dim(mat), g$dim, info = info)
  expect_identical(sum(is.na(mat)), g$nNA, info = info)

  new <- as.vector(mat)[g$sample_idx]
  old <- g$sample_val

  # NA positions must match, and so must the NA-vs-NaN split: a relative
  # comparison alone would treat them as interchangeable.
  expect_identical(is.na(new), is.na(old), info = info)
  expect_identical(is.nan(new), is.nan(old), info = info)

  ok <- is.finite(new) & is.finite(old)
  rel <- abs(new[ok] - old[ok]) / (abs(old[ok]) + 1)
  # expect_lte/expect_lt route `...` to compare(), not to an `info` argument, so
  # use expect_true to keep the magnitude in the failure report.
  max_rel <- max(rel)
  expect_true(max_rel <= GOLDEN_TOL,
              info = sprintf("%s: max rel dev %.3e > %.0e", info, max_rel, GOLDEN_TOL))
  expect_equal(sum(mat, na.rm = TRUE), g$sum, tolerance = GOLDEN_TOL)
}

.run_case <- function(args) {
  fr <- .frames()
  suppressMessages(do.call(
    propaneWarp,
    c(list(image_in = fr$src, keyvalues_out = fr$ref$keyvalues,
           dim_out = fr$dim), args)
  ))
}

CASE_ARGS <- list(
  default      = list(),
  forward      = list(direction = "forward"),
  nearest      = list(interpolation = "nearest"),
  linear       = list(interpolation = "linear"),
  neumann      = list(boundary = "neumann"),
  periodic     = list(boundary = "periodic"),
  nofinenorm   = list(dofinenorm = FALSE),
  notightcrop  = list(dotightcrop = FALSE),
  keepcrop     = list(keepcrop = TRUE),
  extratight   = list(keepcrop = TRUE, extratight = TRUE),
  noscale      = list(doscale = FALSE),
  magzero      = list(magzero_out = 25, magzero_in = 30)
)

test_that("golden fixtures are present and complete", {
  expect_true(file.exists(fixture_path))
  g <- .golden()
  expect_setequal(
    names(g$cases),
    c(names(CASE_ARGS), "warpfield_in", "degenerate_pos", "degenerate_neg",
      "degenerate_reference_field_px")
  )
  for (cs in names(CASE_ARGS)) {
    expect_true(g$cases[[cs]]$status %in% c("ok", "nondeterministic"), info = cs)
  }
})

test_that("1.10.1 forward warp is not reproducible (documented pre-existing)", {
  # Recorded so this is visibly a known issue rather than a silent hole in the
  # golden set. CImg's forward warp scatters into the output without atomics,
  # so repeat runs differ and the result depends on OMP_NUM_THREADS.
  # If this ever starts reporting reproducible, the fixture should be
  # regenerated and the case tightened to bit-identity.
  g <- .golden()$cases$forward
  expect_false(isTRUE(g$reproducible))
})

for (case in names(CASE_ARGS)) {
  local({
    cs <- case
    args <- CASE_ARGS[[cs]]
    test_that(paste0("propaneWarp output matches 1.10.1 within tolerance: ", cs), {
      g <- .golden()$cases[[cs]]
      r <- .run_case(args)

      if (isTRUE(g$reproducible)) {
        .expect_close_golden(r$imDat, g, cs)

        # The sample above is 4096 of ~3.2M pixels. What makes the rest of the
        # image trustworthy is that the run is *stable*: identical digests over
        # repeated runs means the full pixel vector is pinned, so a change
        # anywhere shows up here rather than hiding between sample points.
        d2 <- .digest_matrix(.run_case(args)$imDat)
        expect_identical(d2, .digest_matrix(r$imDat),
                         info = paste(cs, "run-to-run"))
        # Also record against the fixture for human diagnosis; this is *not*
        # asserted, since it is exactly the quantity that legitimately moves
        # when the build environment changes.
        cat(sprintf("[%s] digest %s (1.10.1 was %s)\n", cs, d2, g$digest))
      } else {
        # 1.10.1 itself is not reproducible for this case (forward warp race),
        # so no exact assertion is possible. Fall back to loose checks that the
        # computation has not gone badly wrong. Measured 1.10.1 run-to-run
        # spread: total flux agrees to ~3e-6, and ~0.1% of pixels differ (the
        # race is in CImg's non-atomic forward scatter).
        expect_lt(abs(sum(is.na(r$imDat)) - g$nNA) / g$nNA, g$rel_tol)
        expect_equal(sum(r$imDat, na.rm = TRUE), g$sum, tolerance = g$rel_tol)
        new <- as.vector(r$imDat)[g$sample_idx]
        ok <- is.finite(new) & is.finite(g$sample_val)
        rel <- abs(new[ok] - g$sample_val[ok]) / (abs(g$sample_val[ok]) + 1)
        expect_gte(mean(rel < 1e-6), 0.95)
      }

      kvn <- intersect(names(g$kv), names(r$keyvalues))
      expect_identical(unname(r$keyvalues[kvn]), unname(g$kv[kvn]))
      expect_identical(names(r$keyvalues)[kvn], names(g$kv)[kvn])
    })
  })
}

test_that("a supplied warpfield reproduces the internally built one", {
  fr <- .frames()
  built <- suppressMessages(
    propaneWarp(fr$src, keyvalues_out = fr$ref$keyvalues, dim_out = fr$dim,
                warpfield_return = TRUE)
  )
  wf <- built$warpfield
  expect_s3_class(wf, "cimg")
  r <- .run_case(list(warpfield = wf))

  # The load-bearing comparison is against the run that built the field in the
  # same call: handing back the same field must return the same bits.
  expect_identical(.digest_matrix(r$imDat), .digest_matrix(built$imDat),
                   info = "supplied field != internally built field")

  # The golden record is a second, environment-dependent check.
  g <- .golden()$cases$warpfield_in
  .expect_close_golden(r$imDat, g, "warpfield_in")
})

test_that("a non-overlapping target warp returns blank without building a field", {
  # An inverted tight crop used to be handed straight to the *box* crop, which
  # silently EXPANDS the working grid. At CRPIX1 - 10000 the warp field grew to
  # 8183x1768 = 14.5M px (9.7s, mostly wcslib non-convergence) versus the
  # 3.2M-px reference, and the returned image was a normal-sized block of NaN
  # (0/0 inside the dofinenorm divide) with no explanation. There is no data to
  # warp when the frames do not overlap, so we now skip the crop and the warp
  # and return the requested frame filled with `blank`.
  for (cs in c("degenerate_pos", "degenerate_neg")) {
    g <- .golden()$cases[[cs]]
    fr <- .frames()
    kv_bad <- fr$ref$keyvalues
    kv_bad$CRPIX1 <- kv_bad$CRPIX1 + g$shift

    msgs <- NULL
    t0 <- proc.time()[3]
    # Muffle via the handler itself; wrapping in suppressMessages would consume
    # the condition before this handler ever sees it.
    r <- withCallingHandlers(
      propaneWarp(fr$src, keyvalues_out = kv_bad,
                  dim_out = fr$dim, blank = -99),
      message = function(c) { msgs <<- c(msgs, conditionMessage(c))
                              invokeRestart("muffleMessage") }
    )
    dt <- proc.time()[3] - t0

    # The announced fast path, not a silent one.
    expect_true(any(grepl("does not overlap", msgs)),
                info = paste(cs, paste(msgs, collapse = " | ")))

    # Correct full-size output frame, entirely blank.
    expect_identical(dim(r$imDat), g$out_dim, info = cs)
    expect_identical(sum(is.na(r$imDat)), 0L, info = cs)
    expect_true(all(r$imDat == -99), info = cs)

    # The header describes the frame that was asked for: no expansion and no
    # CRPIX shift, since nothing was cropped.
    kvn <- intersect(names(g$kv), names(r$keyvalues))
    want <- g$kv[kvn]
    want$CRPIX1 <- kv_bad$CRPIX1
    want$XCUTLO <- 1L
    want$XCUTHI <- fr$dim[1]
    want$YCUTLO <- 1L
    want$YCUTHI <- fr$dim[2]
    expect_identical(unname(r$keyvalues[kvn]), unname(want), info = cs)

    # No field was built at all, and the whole call is far cheaper than before.
    expect_null(r$warpfield, info = cs)
    expect_true(dt < 2, info = paste(cs, sprintf("%.2fs", dt)))
  }

  # `blank` is honoured, which the old NaN result never was.
  fr <- .frames()
  kv_bad <- fr$ref$keyvalues
  kv_bad$CRPIX1 <- kv_bad$CRPIX1 - 10000L
  r0 <- suppressMessages(propaneWarp(fr$src, keyvalues_out = kv_bad,
                                     dim_out = fr$dim, blank = 0))
  expect_identical(sum(r0$imDat), 0)
  expect_identical(sum(is.na(r0$imDat)), 0L)
  rNA <- suppressMessages(propaneWarp(fr$src, keyvalues_out = kv_bad,
                                      dim_out = fr$dim))
  expect_identical(as.numeric(sum(is.na(rNA$imDat))), as.numeric(prod(fr$dim)))

  # A warp with genuine overlap must stay quiet and un-shortcut.
  expect_message(suppressMessages(
    propaneWarp(fr$src, keyvalues_out = fr$ref$keyvalues, dim_out = fr$dim)
  ), NA)

  # ...and shifts just inside the overlap boundary still go the normal route.
  for (sh in c(1400L, 1800L)) {
    kvb <- fr$ref$keyvalues
    kvb$CRPIX1 <- kvb$CRPIX1 + sh
    expect_no_warning(suppressMessages(
      r <- propaneWarp(fr$src, keyvalues_out = kvb, dim_out = fr$dim,
                       warpfield_return = TRUE)
    ))
    expect_false(is.null(r$warpfield), info = paste0("sh=", sh))
  }
})

test_that("propaneWarpProPane shares one warp field across bands, bit-identically", {
  fr <- .frames()
  mk <- function() {
    x <- suppressMessages(Rfits_read_image(
      file.path(system.file("extdata/stack", package = "ProPane"),
                .golden()$src_frame)))
    structure(list(image = x, weight = x, inVar = x, exp = x, cold = x),
              class = "ProPane")
  }

  a <- suppressMessages(propaneWarpProPane(mk(), keyvalues_out = fr$ref$keyvalues,
                                           dim_out = fr$dim, warpfield_share = FALSE))
  b <- suppressMessages(propaneWarpProPane(mk(), keyvalues_out = fr$ref$keyvalues,
                                           dim_out = fr$dim, warpfield_share = TRUE))

  expect_named(b, names(a))
  for (n in names(a)) {
    if (is.null(a[[n]])) { expect_null(b[[n]]); next }
    expect_identical(dim(a[[n]]$imDat), dim(b[[n]]$imDat), info = n)
    # The load-bearing assertion: reusing the field must not change a single bit.
    expect_identical(.digest_matrix(a[[n]]$imDat), .digest_matrix(b[[n]]$imDat),
                     info = n)
    expect_identical(a[[n]]$keyvalues, b[[n]]$keyvalues, info = n)
  }

  # Sharing is on by default.
  d <- suppressMessages(propaneWarpProPane(mk(), keyvalues_out = fr$ref$keyvalues,
                                           dim_out = fr$dim))
  expect_identical(.digest_matrix(d$image$imDat), .digest_matrix(b$image$imDat))

  # An explicit warpfield in ... must be honoured, not overridden by sharing.
  wf <- suppressMessages(
    propaneWarp(fr$src, keyvalues_out = fr$ref$keyvalues, dim_out = fr$dim,
                warpfield_return = TRUE))$warpfield
  e <- suppressMessages(propaneWarpProPane(mk(), keyvalues_out = fr$ref$keyvalues,
                                           dim_out = fr$dim, warpfield = wf))
  expect_identical(.digest_matrix(e$image$imDat), .digest_matrix(b$image$imDat))

  # Geometry mismatch forces a rebuild rather than a wrong reuse: give the
  # weight band a different frame's geometry and confirm it still matches what
  # an unshared run produces for that same mixed input.
  #
  # direction="backward" is forced because image_3 has a coarser plate scale,
  # so "auto" would pick the forward warp -- which 1.10.1 cannot reproduce run
  # to run (CImg's non-atomic forward scatter), making any digest comparison
  # between two runs meaningless.
  other <- suppressMessages(Rfits_read_image(
    file.path(system.file("extdata/stack", package = "ProPane"), "image_3.fits")))
  expect_false(identical(dim(other)[1:2], dim(fr$src)[1:2]))

  mixed_shared <- mk(); mixed_shared$weight <- other
  mixed_plain  <- mk(); mixed_plain$weight  <- other

  f <- suppressMessages(propaneWarpProPane(mixed_shared, keyvalues_out = fr$ref$keyvalues,
                                           dim_out = fr$dim, warpfield_share = TRUE,
                                           direction = "backward"))
  h <- suppressMessages(propaneWarpProPane(mixed_plain, keyvalues_out = fr$ref$keyvalues,
                                           dim_out = fr$dim, warpfield_share = FALSE,
                                           direction = "backward"))

  # The mismatched band must be rebuilt, not warped with the wrong field...
  expect_identical(.digest_matrix(f$weight$imDat), .digest_matrix(h$weight$imDat))
  # ...and every other band must still be bit-identical.
  for (n in setdiff(names(f), "weight")) {
    if (is.null(f[[n]])) next
    expect_identical(.digest_matrix(f[[n]]$imDat), .digest_matrix(h[[n]]$imDat),
                     info = n)
  }
})

test_that("warpgrid defaults to the exact path", {
  expect_true("warpgrid" %in% names(formals(propaneWarp)))
  expect_identical(formals(propaneWarp)$warpgrid, "exact")
  expect_identical(formals(propaneWarp)$warptol, 1e-05)
  # Default behaviour remains within tolerance of the golden record.
  g <- .golden()$cases$default
  .expect_close_golden(.run_case(list())$imDat, g, "warpgrid=default")
})

test_that("coarse warpgrid is not the same as the exact path", {
  # Opt-in and explicitly approximate, so it is *expected* to differ from the
  # golden record -- but only by the small amount a verified reconstruction
  # should, never by a pixel or more.
  skip_if_not("warpgrid" %in% names(formals(propaneWarp)))
  g <- .golden()$cases$default
  fr <- .frames()

  r <- suppressMessages(propaneWarp(fr$src, keyvalues_out = fr$ref$keyvalues,
                                    dim_out = fr$dim, warpgrid = "approx",
                                    warpfield_return = TRUE))
  expect_s3_class(r$warpfield, "cimg")

  # The field is a *different object* from the exact one, but must agree to far
  # better than a pixel. Compare against the field the exact path produces for
  # the same tight-cropped geometry.
  ex <- suppressMessages(propaneWarp(fr$src, keyvalues_out = fr$ref$keyvalues,
                                     dim_out = fr$dim, warpfield_return = TRUE))
  a <- cbind(as.vector(ex$warpfield[, , 1]), as.vector(ex$warpfield[, , 2]))
  b <- cbind(as.vector(r$warpfield[, , 1]), as.vector(r$warpfield[, , 2]))
  expect_identical(dim(a), dim(b))
  ok <- is.finite(a) & is.finite(b)
  expect_lt(max(abs(a[ok] - b[ok])), 1e-3)

  # Not bit-identical to exact (it is an approximation), but must agree
  # photometrically: integrated flux to ~1e-5, NA pattern unchanged.
  expect_identical(dim(r$imDat), g$dim)
  expect_equal(sum(r$imDat, na.rm = TRUE), g$sum, tolerance = 1e-5)
  expect_identical(sum(is.na(r$imDat)), g$nNA)
})

test_that("coarse warpgrid refuses rather than returning an inaccurate field", {
  skip_if_not(exists(".warpfield_coarse", envir = asNamespace("ProPane")))
  coarse <- get(".warpfield_coarse", envir = asNamespace("ProPane"))
  out2in <- get(".warpfunc_out2in", envir = asNamespace("ProPane"))

  suppressMessages(library(Rfits))
  src <- suppressMessages(Rfits_read_image(
    file.path(system.file("extdata/stack", package = "ProPane"), "image_1.fits")))
  ref <- suppressMessages(Rfits_read_image(
    file.path(system.file("extdata/stack", package = "ProPane"), "image_2.fits")))
  hi <- Rfits_keyvalues_to_raw(src$keyvalues)
  ho <- Rfits_header_to_raw(Rfits_keyvalues_to_header(ref$keyvalues))
  wf <- function(x, y, cores = 1, ...)
    out2in(x, y, header_in = hi, header_out = ho, cores = cores)

  # An absurdly tight tolerance cannot be met at step 2, so the builder must
  # return NULL and let the caller fall back to the exact path.
  expect_null(coarse(wf, c(400, 400), tol = 1e-30, step0 = 64L, max_refine = 2L))

  # A grid too small to interpolate must also refuse.
  expect_null(coarse(wf, c(3, 3), tol = 1e-5, step0 = 64L))

  # And a reachable tolerance must succeed, with the claimed error agreeing
  # with the true error against the full-resolution field.
  r <- suppressMessages(coarse(wf, c(400, 400), tol = 1e-5, step0 = 64L))
  expect_false(is.null(r))
  pg <- expand.grid(seq_len(400), seq_len(400))
  ref_f <- suppressWarnings(wf(pg[, 1], pg[, 2], cores = 1))
  got <- cbind(as.vector(r$warpfield[, , 1]), as.vector(r$warpfield[, , 2]))
  ok <- is.finite(ref_f) & is.finite(got)
  true_err <- max(abs(got[ok] - ref_f[ok]))
  expect_lt(true_err, 1e-3)
  # The self-reported midpoint error should be the same order as the truth:
  # it is an estimate, so allow a generous factor either way.
  expect_gt(r$maxerr, true_err / 100)
  expect_lt(r$maxerr, true_err * 100)
})
