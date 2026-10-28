# Golden-fixture regression tests for propaneWarp().
#
# fixtures/warp-golden.rds is regenerated with bench/make-fixtures.R, which must
# be run against an installed build of the code whose behaviour is the target.
# It was originally produced by released ProPane 1.10.1; the current record was
# produced by 1.12.1 so that it stores the ZIMAGE size keys (see below) and the
# deterministic forward warp. Steps 1-3 of the warp performance work must not
# change the output; only the opt-in `warpgrid` path is permitted to differ.
#
# "Must not change" is enforced to a tight tolerance rather than bit-for-bit,
# because exact digests do not survive a change of build environment. The
# original 1.10.1 record showed a *uniform* residual of ~8.4e-14 against a build
# in a different Rfits/Rwcs environment, traceable to pixscale(): its CD matrix
# round-trips through 15-significant-digit FITS cards, and it is sampled at the
# field centre, where a 0.001 px shift already moves the scale by 2e-9.
#
# Now that the record and the code come from the same build, the sampled pixels
# agree exactly on the machine that generated it (max relative deviation 0), and
# the printed digests match the fixture. The tolerance is kept anyway: it is what
# lets a contributor on a different wcslib build get a useful failure, and 1e-9
# is still far below the O(1e-3) that a real warp or crop bug moves pixels by.
#
# So the gates below are: exact shape, exact missing-value pattern (NA count, NA
# positions, NA vs NaN), every sampled pixel agreeing to GOLDEN_TOL relative, and
# bit-stability across repeated runs.
#
# COST. The full suite used to take ~110s here and roughly 4-5x that on the
# GitHub Actions Ubuntu runners, which pushed the whole check past the point
# where it was worth waiting for. Two changes bring it down to ~35s without
# giving up any case:
#
#   - the run-to-run bit-stability re-run is now done for `default` only. It
#     proves the warp is deterministic, which is a property of the code path
#     rather than of each option, so paying for it twelve times bought nothing.
#     The fixture still records every case's digest for human diagnosis.
#   - `extratight` is no longer run separately: at this geometry it produces a
#     byte-identical result to `keepcrop` (same digest in the fixture), so it
#     cannot fail without `keepcrop` failing first.
#
# The `propaneWarpProPane` warpfield-sharing test was the single most expensive
# thing in the file (~46s for a dozen full-size warps). It is kept, but at a
# smaller output frame -- its assertions compare two runs of the same call to
# each other and need no golden record, so shrinking dim_out costs no coverage.

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

# Compare a result's header against the golden record. The fixture stores both
# NAXIS* and ZNAXIS*, so this is a straight comparison: on a ZIMAGE frame the
# image size lives in ZNAXIS1/2 and NAXIS1/2 are the compressed BINTABLE row
# width (bytes) and row count, and having both pairs recorded is what makes that
# distinction visible rather than worked around.
.expect_kv_golden <- function(r, g, info) {
  kvn <- intersect(names(g$kv), names(r$keyvalues))
  expect_identical(unname(r$keyvalues[kvn]), unname(g$kv[kvn]), info = info)
  expect_identical(names(r$keyvalues)[kvn], names(g$kv)[kvn], info = info)

  # Whichever pair is authoritative for this frame must describe the matrix.
  if (isTRUE(r$keyvalues$ZIMAGE)) {
    expect_identical(c(r$keyvalues$ZNAXIS1, r$keyvalues$ZNAXIS2),
                     as.integer(dim(r$imDat)), info = paste(info, "ZNAXIS vs dim"))
  } else {
    expect_identical(c(r$keyvalues$NAXIS1, r$keyvalues$NAXIS2),
                     as.integer(dim(r$imDat)), info = paste(info, "NAXIS vs dim"))
  }
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
  noscale      = list(doscale = FALSE),
  magzero      = list(magzero_out = 25, magzero_in = 30)
)

test_that("golden fixtures are present and complete", {
  expect_true(file.exists(fixture_path))
  g <- .golden()
  # `extratight` is in the fixture but not in CASE_ARGS: it is byte-identical to
  # `keepcrop` at this geometry, so running it separately only doubled the cost
  # of a case that cannot fail on its own.
  expect_setequal(
    names(g$cases),
    c(names(CASE_ARGS), "extratight", "warpfield_in", "degenerate_pos",
      "degenerate_neg", "degenerate_reference_field_px")
  )
  for (cs in names(CASE_ARGS)) {
    expect_true(g$cases[[cs]]$status %in% c("ok", "nondeterministic"), info = cs)
  }
})

test_that("the forward warp is reproducible at this geometry", {
  # Under 1.10.1 the forward warp was recorded as NON-reproducible: CImg scatters
  # into the output without atomics, so repeat runs differed and the result was
  # said to depend on OMP_NUM_THREADS. That is why the fixture used to carry
  # status "nondeterministic" for this case and why the assertions fell back to
  # loose flux checks instead of bit comparison.
  #
  # On the current build it is deterministic: 6 consecutive runs give one digest,
  # and OMP_NUM_THREADS of 1, 2, 4 and 8 all give the same one. So the case is now
  # recorded as reproducible and checked exactly like the others. If it ever
  # flakes again, this test is the signal -- it will fail rather than quietly
  # resting on the looser fallback.
  g <- .golden()$cases$forward
  expect_true(isTRUE(g$reproducible))

  # Independence from the thread count is the property that made 1.10.1 fail, so
  # it is worth asserting directly rather than trusting the fixture record.
  fr <- .frames()
  old <- Sys.getenv("OMP_NUM_THREADS", unset = NA)
  on.exit(if (is.na(old)) Sys.unsetenv("OMP_NUM_THREADS") else Sys.setenv(OMP_NUM_THREADS = old))
  digs <- vapply(c("1", "2", "4"), function(n) {
    Sys.setenv(OMP_NUM_THREADS = n)
    .digest_matrix(.run_case(list(direction = "forward"))$imDat)
  }, character(1))
  expect_identical(length(unique(digs)), 1L,
                   info = paste(names(digs), digs, collapse = " | "))
})

for (case in names(CASE_ARGS)) {
  local({
    cs <- case
    args <- CASE_ARGS[[cs]]
    test_that(paste0("propaneWarp output matches the golden record: ", cs), {
      g <- .golden()$cases[[cs]]
      r <- .run_case(args)

      if (isTRUE(g$reproducible)) {
        .expect_close_golden(r$imDat, g, cs)

        # The sample above is 4096 of ~3.2M pixels. What makes the rest of the
        # image trustworthy is that the run is *stable*: identical digests over
        # repeated runs means the full pixel vector is pinned, so a change
        # anywhere shows up here rather than hiding between sample points.
        #
        # This costs a second full warp, so it is done for `default` only.
        # Determinism is a property of the warp code path, not of whichever
        # interpolation or boundary option happens to be set, so proving it once
        # is enough -- running it for all eleven cases roughly doubled the file
        # without widening what a regression could be caught by.
        if (cs == "default") {
          d2 <- .digest_matrix(.run_case(args)$imDat)
          expect_identical(d2, .digest_matrix(r$imDat),
                           info = paste(cs, "run-to-run"))
        }
        # Also record against the fixture for human diagnosis; this is *not*
        # asserted, since it is exactly the quantity that legitimately moves
        # when the build environment changes.
        cat(sprintf("[%s] digest %s (fixture %s)\n",
                    cs, .digest_matrix(r$imDat), g$digest))
      } else {
        # A case recorded as non-reproducible cannot be asserted exactly, so
        # fall back to loose checks that the computation has not gone badly
        # wrong: total flux agrees to ~3e-6 and ~0.1% of pixels differ. No case
        # is in this state on the current fixture (the forward warp is now
        # deterministic, see above); the branch is kept so a future
        # nondeterministic case is handled deliberately rather than silently
        # compared bit-for-bit and left flaking.
        expect_lt(abs(sum(is.na(r$imDat)) - g$nNA) / g$nNA, g$rel_tol)
        expect_equal(sum(r$imDat, na.rm = TRUE), g$sum, tolerance = g$rel_tol)
        new <- as.vector(r$imDat)[g$sample_idx]
        ok <- is.finite(new) & is.finite(g$sample_val)
        rel <- abs(new[ok] - g$sample_val[ok]) / (abs(g$sample_val[ok]) + 1)
        expect_gte(mean(rel < 1e-6), 0.95)
      }

      .expect_kv_golden(r, g, cs)
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
    # NAXIS1/2 are the compressed-storage geometry on these frames and are not
    # touched by the crop, so they come straight from the fixture; the frame that
    # was asked for is the full input, recorded in ZNAXIS1/2.
    expect_identical(unname(r$keyvalues[kvn]), unname(want), info = cs)
    expect_identical(c(r$keyvalues$ZNAXIS1, r$keyvalues$ZNAXIS2),
                     as.integer(fr$dim), info = cs)

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
  # Every assertion below compares two runs of the *same* call to each other, so
  # none of it depends on the golden record and none needs the full 1816x1767
  # frame. At 700x700 the warp still covers real sky (99.7% finite, 1400 NA)
  # and the whole test costs about half as much.
  dm <- c(700L, 700L)
  mk <- function() {
    x <- suppressMessages(Rfits_read_image(
      file.path(system.file("extdata/stack", package = "ProPane"),
                .golden()$src_frame)))
    structure(list(image = x, weight = x, inVar = x, exp = x, cold = x),
              class = "ProPane")
  }

  a <- suppressMessages(propaneWarpProPane(mk(), keyvalues_out = fr$ref$keyvalues,
                                           dim_out = dm, warpfield_share = FALSE))
  b <- suppressMessages(propaneWarpProPane(mk(), keyvalues_out = fr$ref$keyvalues,
                                           dim_out = dm, warpfield_share = TRUE))

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
                                           dim_out = dm))
  expect_identical(.digest_matrix(d$image$imDat), .digest_matrix(b$image$imDat))

  # An explicit warpfield in ... must be honoured, not overridden by sharing.
  wf <- suppressMessages(
    propaneWarp(fr$src, keyvalues_out = fr$ref$keyvalues, dim_out = dm,
                warpfield_return = TRUE))$warpfield
  e <- suppressMessages(propaneWarpProPane(mk(), keyvalues_out = fr$ref$keyvalues,
                                           dim_out = dm, warpfield = wf))
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
                                           dim_out = dm, warpfield_share = TRUE,
                                           direction = "backward"))
  h <- suppressMessages(propaneWarpProPane(mixed_plain, keyvalues_out = fr$ref$keyvalues,
                                           dim_out = dm, warpfield_share = FALSE,
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
  expect_identical(formals(propaneWarp)$warptol, 1e-03)
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

# ---------------------------------------------------------------------------
# Degenerate output frames: 1x1, 1xN and Nx1.
#
# imager::warp() returns a `cimg`, and as.matrix() does not preserve the axes of
# a cimg that is one pixel wide. as.matrix.cimg() has two relevant branches:
#
#   - exactly two dimensions equal to 1 (e.g. 1 x N) -> squeeze() then
#     class<- "matrix", giving an N x 1 matrix;
#   - exactly one dimension greater than 1 (e.g. 1 x N or N x 1) -> warns
#     "Image is one-dimensional" and returns as.vector(x) %>% as.matrix, which
#     is ALWAYS an N x 1 matrix;
#   - a genuinely 1 x 1 x 1 x 1 cimg -> error "Too many non-empty dimensions".
#
# So a 1 x N frame could come back as N x 1, and a 1 x 1 frame could not be
# converted at all. The conversion after the warp is therefore an explicit
#   matrix(x, nrow = dim(x)[1], ncol = dim(x)[2])
# which rebuilds the axes the cimg carried instead of inferring them.
#
# WHY keepcrop / dotightcrop MATTER. The closing crop at the end of
# propaneWarp() subsets with box=1, i.e. it rebuilds the matrix from the crop
# range rather than from the warped matrix. That repairs a transposed
# intermediate, so the default path masked the bug and it only reached callers
# through the uncropped routes: keepcrop = TRUE, or dotightcrop = FALSE. These
# tests pin both routes; the uncropped ones are the ones that can actually fail.

# Shared helpers for the degenerate-shape tests. A one-pixel-wide input plus an
# integer CRPIX shift is what drags the working grid down to a single column (or
# row), which is the condition that makes as.matrix.cimg() transpose. Without the
# shift the input is padded and the frame never becomes one pixel wide, so a
# future change could silently stop exercising the bug.
.degenerate_warp <- function(xr, yr, dim_out, shift = c(0L, 0L), ...) {
  fr <- .frames()
  im <- suppressMessages(fr$src[xr, yr])
  im$imDat[is.na(im$imDat)] <- 0
  kv <- im$keyvalues
  kv$NAXIS1 <- dim_out[1]
  kv$NAXIS2 <- dim_out[2]
  kv$CRPIX1 <- kv$CRPIX1 + shift[1]
  kv$CRPIX2 <- kv$CRPIX2 + shift[2]
  suppressMessages(propaneWarp(
    im, keyvalues_out = kv, dim_out = dim_out, direction = "backward",
    warpfield_return = TRUE, ...
  ))
}

# The frame is returned at `want`, with axes in that order. For the uncropped
# routes `want` is the request; the default cropped route returns the crop, which
# is asserted separately.
.expect_shape <- function(r, want, info) {
  # Two dimensions, not a bare vector.
  expect_identical(length(dim(r$imDat)), 2L, info = info)
  # This is the assertion the old as.matrix() conversion failed: 1 x N came back
  # as N x 1.
  expect_identical(dim(r$imDat), as.integer(want), info = info)
  # The FITS header must describe the same frame as the matrix, so a transposed
  # imDat is caught even by a caller who never looks at dim(). XCUT* are numeric
  # here and integer elsewhere, hence expect_equal.
  expect_equal(
    c(r$keyvalues$XCUTHI - r$keyvalues$XCUTLO + 1,
      r$keyvalues$YCUTHI - r$keyvalues$YCUTLO + 1),
    as.numeric(dim(r$imDat)),
    info = info
  )
}

test_that("propaneWarp returns a 1x1 frame rather than failing to convert it", {
  # imager's 1 x 1 x 1 x 1 branch: as.matrix.cimg() used to error out here with
  # "Too many non-empty dimensions" rather than return anything at all.
  r <- .degenerate_warp(c(900L, 900L), c(900L, 900L), c(1L, 1L),
                        dotightcrop = FALSE)
  .expect_shape(r, c(1L, 1L), "1x1 dotightcrop=FALSE")
  expect_identical(dim(r$warpfield)[1:2], c(1L, 1L))
})

test_that("propaneWarp keeps a 1xN frame 1 wide, not transposed to Nx1", {
  # The reported case. Old conversion returned 20 x 1 while the header said
  # 1 x 20.
  r <- .degenerate_warp(c(900L, 900L), c(900L, 919L), c(1L, 20L),
                        shift = c(-1L, 0L), keepcrop = TRUE)
  .expect_shape(r, c(1L, 20L), "1x20 keepcrop")

  # Same shape via the other uncropped route.
  r <- .degenerate_warp(c(900L, 900L), c(900L, 919L), c(1L, 20L),
                        shift = c(-1L, 0L), keepcrop = TRUE, dotightcrop = FALSE)
  .expect_shape(r, c(1L, 20L), "1x20 keepcrop, dotightcrop=FALSE")
})

test_that("propaneWarp keeps an Nx1 frame 1 tall", {
  # The mirror image. Already correct before the fix, and kept so a "fix" that
  # transposes everything cannot pass.
  r <- .degenerate_warp(c(900L, 909L), c(900L, 900L), c(10L, 1L),
                        shift = c(-1L, -1L), keepcrop = TRUE)
  .expect_shape(r, c(10L, 1L), "10x1 keepcrop")
})

test_that("the default cropped route stays self-consistent at degenerate shapes", {
  # With the closing crop the returned frame is the crop, not the request, and
  # subsetting with box=1 rebuilds the axes -- which is exactly what hid this bug
  # on the default path. Only header/matrix agreement is meaningful here.
  for (dm in list(c(1L, 1L), c(1L, 10L), c(10L, 1L))) {
    r <- .degenerate_warp(c(900L, 900L), c(900L, 909L), dm)
    info <- paste("default", paste(dm, collapse = "x"))
    expect_identical(length(dim(r$imDat)), 2L, info = info)
    expect_equal(
      c(r$keyvalues$XCUTHI - r$keyvalues$XCUTLO + 1,
        r$keyvalues$YCUTHI - r$keyvalues$YCUTLO + 1),
      as.numeric(dim(r$imDat)),
      info = info
    )
  }
})

test_that("a 1-pixel-wide warp keeps the 1xN shape for every interpolation", {
  # The transposition lived in the cimg -> matrix conversion, which every
  # interpolation and boundary mode shares, so one geometry per option is enough
  # to show none of them reaches the return value by a different route.
  for (interp in c("nearest", "linear", "cubic")) {
    for (bound in c("dirichlet", "neumann", "periodic")) {
      r <- .degenerate_warp(c(900L, 900L), c(900L, 919L), c(1L, 20L),
                            shift = c(-1L, 0L), keepcrop = TRUE,
                            interpolation = interp, boundary = bound)
      info <- paste(interp, bound)
      expect_identical(dim(r$imDat), c(1L, 20L), info = info)
      expect_identical(length(dim(r$imDat)), 2L, info = info)
    }
  }
})

# ---------------------------------------------------------------------------
# propaneRebin carries the same cimg -> matrix conversion, and had the same
# defect: as.matrix() turned a 1 x N rebin into N x 1 and errored on 1 x 1.
# Unlike the warp there is a second trap here -- the output grid is the RESIZED
# one, not the input one. An explicit
#   matrix(x, nrow = dim(input)[1], ncol = dim(input)[2])
# would compile happily and silently undo every scale != 1, so the shape tests
# below are paired with assertions that the resized grid is still the resized
# grid.

test_that("propaneRebin keeps the resized grid, not the input grid", {
  fr <- .frames()
  im <- suppressMessages(fr$src[c(900L, 929L), c(900L, 919L)])  # 30 x 20
  im$imDat[is.na(im$imDat)] <- 0
  expect_identical(dim(im$imDat), c(30L, 20L))

  # Same size_x/size_y arithmetic the function uses.
  for (sc in c(0.5, 2, 3)) {
    want <- if (sc > 1) {
      floor(c(30, 20) * sc - (sc - 1))
    } else {
      floor(c(30, 20) * sc)
    }
    r <- suppressWarnings(propaneRebin(im, scale = sc))
    expect_identical(dim(r$imDat), as.integer(want), info = paste("scale", sc))
    # These frames are ZIMAGE, so the image size lives in ZNAXIS*.
    expect_identical(c(r$keyvalues$ZNAXIS1, r$keyvalues$ZNAXIS2), as.integer(want),
                     info = paste("scale", sc))
  }
})

test_that("propaneRebin leaves NAXIS1/2 alone on a ZIMAGE frame", {
  # On a tile-compressed image NAXIS1/2 are NOT the image size. In the FITS
  # convention used here (RICE_1, 8-byte doubles) the extension is a BINTABLE
  # whose NAXIS1 is the compressed row width in BYTES and whose NAXIS2 is the
  # table row count; the real image lives in ZNAXIS1/2. For the packaged
  # image_1.fits that is NAXIS1 = 24 (tile payload), NAXIS2 = 1836 (rows), while
  # the image is 1816 x 1836.
  #
  # So "keeping the header in step" on a ZIMAGE frame means rewriting ZNAXIS*
  # and leaving NAXIS* untouched: writing the image size into NAXIS1/2 would
  # describe a table that is nowhere on disk. This is the behaviour of
  # Rfits_check_image(), of propaneWarp(), and of propaneRebin() -- the test
  # pins all three together so they cannot drift apart silently.
  fr <- .frames()
  im <- suppressMessages(fr$src[c(900L, 929L), c(900L, 919L)])  # 30 x 20
  im$imDat[is.na(im$imDat)] <- 0
  expect_true(isTRUE(im$keyvalues$ZIMAGE))

  # The storage keys on the input are not the image size, which is the whole
  # point: this fixture would fail if they happened to coincide.
  expect_false(identical(c(im$keyvalues$NAXIS1, im$keyvalues$NAXIS2), c(30L, 20L)))

  for (sc in c(0.5, 2)) {
    info <- paste("scale", sc)
    want <- if (sc > 1) floor(c(30, 20) * sc - (sc - 1)) else floor(c(30, 20) * sc)

    r <- suppressWarnings(propaneRebin(im, scale = sc))
    expect_identical(dim(r$imDat), as.integer(want), info = info)
    # Image size goes to ZNAXIS*...
    expect_identical(c(r$keyvalues$ZNAXIS1, r$keyvalues$ZNAXIS2), as.integer(want),
                     info = info)
    # ...and the compressed-storage keys are left exactly as they were.
    expect_identical(c(r$keyvalues$NAXIS1, r$keyvalues$NAXIS2),
                     c(im$keyvalues$NAXIS1, im$keyvalues$NAXIS2),
                     info = info)
  }

  # propaneWarp follows the same rule, whether dim_out is derived from the header
  # or supplied (the supplied-dim_out case is covered in its own test below).
  w <- suppressMessages(propaneWarp(im, keyvalues_out = im$keyvalues,
                                    direction = "backward"))
  expect_identical(dim(w$imDat), c(30L, 20L))
  expect_identical(c(w$keyvalues$ZNAXIS1, w$keyvalues$ZNAXIS2), c(30L, 20L))
  expect_identical(c(w$keyvalues$NAXIS1, w$keyvalues$NAXIS2),
                   c(im$keyvalues$NAXIS1, im$keyvalues$NAXIS2))

  # The non-ZIMAGE branch has no ZNAXIS* at all, so there NAXIS* IS the image
  # size and must be updated. Strip the compression keys to get such a frame.
  plain <- im
  plain$keyvalues$ZIMAGE <- NULL
  for (k in c("ZNAXIS", "ZNAXIS1", "ZNAXIS2", "ZBITPIX", "ZTILE1", "ZTILE2",
              "ZCMPTYPE", "ZNAME1", "ZNAME2", "ZVAL1", "ZVAL2")) {
    plain$keyvalues[[k]] <- NULL
  }
  plain$keynames <- names(plain$keyvalues)
  plain <- Rfits::Rfits_check_image(plain)
  expect_false(isTRUE(plain$keyvalues$ZIMAGE))

  rp <- suppressWarnings(propaneRebin(plain, scale = 0.5))
  expect_identical(dim(rp$imDat), c(15L, 10L))
  expect_identical(c(rp$keyvalues$NAXIS1, rp$keyvalues$NAXIS2), c(15L, 10L))
})

test_that("propaneWarp writes a supplied dim_out to ZNAXIS on a ZIMAGE frame", {
  # Counterpart of the propaneRebin test above, for the branch that takes
  # dim_out from the caller rather than from the header. This branch used to
  # write dim_out into NAXIS1/2 unconditionally, which contradicted its own
  # sibling branch a few lines earlier: with dim_out = NULL the function reads
  # the size back out of ZNAXIS*, so the two routes disagreed about where the
  # image size lives. The result was a header claiming a 15 x 10 image alongside
  # a 15 x 10 compressed table whose ZTILE1 was still 1816.
  fr <- .frames()
  im <- suppressMessages(fr$src[c(900L, 929L), c(900L, 919L)])  # 30 x 20
  im$imDat[is.na(im$imDat)] <- 0
  expect_true(isTRUE(im$keyvalues$ZIMAGE))
  storage <- c(im$keyvalues$NAXIS1, im$keyvalues$NAXIS2)
  expect_false(identical(storage, c(30L, 20L)))

  # Supplying dim_out must agree with letting it be derived.
  derived <- suppressMessages(propaneWarp(im, keyvalues_out = im$keyvalues,
                                          direction = "backward"))
  for (dm in list(c(15L, 10L), c(30L, 20L), c(10L, 30L))) {
    info <- paste("dim_out", paste(dm, collapse = "x"))
    w <- suppressMessages(propaneWarp(im, keyvalues_out = im$keyvalues,
                                      dim_out = dm, direction = "backward"))

    expect_identical(dim(w$imDat), dm, info = info)
    # Image size goes to ZNAXIS*...
    expect_identical(c(w$keyvalues$ZNAXIS1, w$keyvalues$ZNAXIS2), dm, info = info)
    # ...and the compressed-storage keys are untouched.
    expect_identical(c(w$keyvalues$NAXIS1, w$keyvalues$NAXIS2), storage, info = info)
  }

  # The two routes must leave the header in the same shape, which is the
  # inconsistency that made the bug visible in the first place.
  w <- suppressMessages(propaneWarp(im, keyvalues_out = im$keyvalues,
                                    dim_out = c(30L, 20L), direction = "backward"))
  expect_identical(c(w$keyvalues$ZNAXIS1, w$keyvalues$ZNAXIS2),
                   c(derived$keyvalues$ZNAXIS1, derived$keyvalues$ZNAXIS2))
  expect_identical(c(w$keyvalues$NAXIS1, w$keyvalues$NAXIS2),
                   c(derived$keyvalues$NAXIS1, derived$keyvalues$NAXIS2))

  # A ZNAXIS-derived round trip still works: what was written can be read back.
  again <- suppressMessages(propaneWarp(im, keyvalues_out = w$keyvalues,
                                        direction = "backward"))
  expect_identical(dim(again$imDat), c(30L, 20L))

  # Non-ZIMAGE frames are unaffected: there NAXIS* is the image size.
  plain <- im
  plain$keyvalues$ZIMAGE <- NULL
  for (k in c("ZNAXIS", "ZNAXIS1", "ZNAXIS2", "ZBITPIX", "ZTILE1", "ZTILE2",
              "ZCMPTYPE", "ZNAME1", "ZNAME2", "ZVAL1", "ZVAL2")) {
    plain$keyvalues[[k]] <- NULL
  }
  plain$keynames <- names(plain$keyvalues)
  plain <- Rfits::Rfits_check_image(plain)

  wp <- suppressMessages(propaneWarp(plain, keyvalues_out = plain$keyvalues,
                                     dim_out = c(15L, 10L), direction = "backward"))
  expect_identical(dim(wp$imDat), c(15L, 10L))
  expect_identical(c(wp$keyvalues$NAXIS1, wp$keyvalues$NAXIS2), c(15L, 10L))
})

test_that("readers of the output header prefer ZNAXIS on a ZIMAGE frame", {
  # Everything in the package that sizes a frame from a header keys off
  # ZNAXIS1/2 when ZIMAGE is set, and falls back to NAXIS1/2 otherwise:
  # propaneWarp's dim_out derivation, propaneStackWarpInVar, propaneWarpDump,
  # propaneStackWarpMed and propaneStackWarpQuanWt. That is the mirror image of
  # the writer rule pinned above -- and the two only agree if both sides use the
  # same pair, which is why the writer bug was visible at all.
  #
  # The fixture makes the two pairs deliberately disagree (image 30 x 20, storage
  # keys 7 x 99), so a reader that silently preferred NAXIS* would come back a
  # wildly different shape rather than off by one pixel.
  fr <- .frames()
  im <- suppressMessages(fr$src[c(900L, 929L), c(900L, 919L)])  # 30 x 20
  im$imDat[is.na(im$imDat)] <- 0
  kv <- im$keyvalues
  kv$ZNAXIS1 <- 30L; kv$ZNAXIS2 <- 20L
  kv$NAXIS1 <- 7L;   kv$NAXIS2 <- 99L
  expect_true(isTRUE(kv$ZIMAGE))

  # Writer: dim_out supplied goes to ZNAXIS*, storage keys untouched.
  w <- suppressMessages(propaneWarp(im, keyvalues_out = kv, dim_out = c(15L, 10L),
                                    direction = "backward"))
  expect_identical(c(w$keyvalues$ZNAXIS1, w$keyvalues$ZNAXIS2), c(15L, 10L))
  expect_identical(c(w$keyvalues$NAXIS1, w$keyvalues$NAXIS2), c(7L, 99L))

  # Reader: with dim_out derived from the header, ZNAXIS wins.
  d <- suppressMessages(propaneWarp(im, keyvalues_out = kv, direction = "backward"))
  expect_identical(dim(d$imDat), c(30L, 20L))

  # Reader: the stacking entry point that sizes its accumulator from the header.
  s <- suppressMessages(propaneStackWarpInVar(
    image_list = list(im), keyvalues_out = kv, dim_out = NULL,
    cores = 1, cores_warp = 1, multitype = "fork"
  ))
  expect_identical(dim(s$image$imDat), c(30L, 20L))

  # Reader: the dump path, which sizes its frames the same way. dim_out = NULL
  # forces it to consult the header rather than the caller.
  out_dir <- file.path(tempdir(), paste0("pp_znaxis_", sample.int(1e6, 1)))
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  suppressMessages(propaneWarpDump(
    image_list = list(im), keyvalues_out = kv, dim_out = NULL,
    cores = 1, cores_warp = 1, dump_dir = out_dir, stub = "z"
  ))
  dumped <- list.files(out_dir, pattern = "\\.fits$", full.names = TRUE)
  expect_gt(length(dumped), 0L)
  back <- suppressMessages(Rfits::Rfits_read_image(dumped[1], ext = 1))
  expect_identical(dim(back$imDat), c(30L, 20L))
})

test_that("propaneRebin keeps 1x1, 1xN and Nx1 shapes", {
  fr <- .frames()
  im <- suppressMessages(fr$src[c(900L, 919L), c(900L, 919L)])  # 20 x 20 source
  im$imDat[is.na(im$imDat)] <- 0

  cases <- list(
    list(name = "1x1",  rows = c(900L, 900L), cols = c(900L, 900L), scale = 1,  want = c(1L, 1L)),
    list(name = "1x10", rows = c(900L, 900L), cols = c(900L, 909L), scale = 1,  want = c(1L, 10L)),
    list(name = "10x1", rows = c(900L, 909L), cols = c(900L, 900L), scale = 1,  want = c(10L, 1L)),
    # scale != 1: a 1-wide input must stay 1 wide and grow along the other axis.
    list(name = "1x10 x2", rows = c(900L, 900L), cols = c(900L, 909L), scale = 2, want = c(1L, 19L))
  )

  for (cs in cases) {
    s <- suppressMessages(fr$src[cs$rows, cs$cols])
    s$imDat[is.na(s$imDat)] <- 0
    r <- suppressWarnings(propaneRebin(s, scale = cs$scale))
    # Exactly two dimensions, in the requested order.
    expect_identical(length(dim(r$imDat)), 2L, info = cs$name)
    expect_identical(dim(r$imDat), cs$want, info = cs$name)
    # Header agrees with the matrix (ZIMAGE frames, so ZNAXIS*).
    expect_identical(c(r$keyvalues$ZNAXIS1, r$keyvalues$ZNAXIS2), cs$want, info = cs$name)
  }

  # The plain-matrix branch shares the conversion and has the same expectations.
  for (dims in list(c(1L, 1L), c(1L, 10L), c(10L, 1L))) {
    m <- matrix(seq_len(prod(dims)), dims[1], dims[2])
    expect_identical(dim(suppressWarnings(propaneRebin(m, scale = 1))), dims,
                     info = paste("matrix", paste(dims, collapse = "x")))
  }
})
