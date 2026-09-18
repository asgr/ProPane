# Accuracy gate for the coarse warp-field path.
#
# The bundled stack frames are a small (~10 arcmin) near-affine field, which
# almost any smooth approximation passes. These synthetic WCS headers stress the
# bilinear reconstruction with the things that could break it: large spherical
# curvature, strong SIP/TPV distortion, high declination, steep rotation, and
# mismatched projections.
#
# The property being tested is NOT "the coarse field is always accurate" -- it is
# the one that actually matters: the builder either achieves its claimed
# tolerance or refuses and lets the caller use the exact field. A silently wrong
# warp would be a bug; a refused warp is only a lost speedup.

test_that("coarse warpgrid is accurate or refuses, across hard WCS geometries", {
  skip_on_cran()
  suppressMessages(library(Rfits))
  suppressMessages(library(Rwcs))
  suppressMessages(library(imager))

  ns <- asNamespace("ProPane")
  coarse <- get(".warpfield_coarse", envir = ns)
  out2in <- get(".warpfunc_out2in", envir = ns)

  # Synthesise a WCS keyvalue list. propaneGenWCS() insists on scanning real
  # files, which is the wrong dependency for a synthetic test.
  mk = function(RA, Dec, pixscale, dim, rot = 0,
                ctype1 = "RA---TAN", ctype2 = "DEC--TAN", extra = list()) {
    ang = pixscale / 3600
    t = rot * pi / 180
    kv = c(
      list(WCSAXES = 2L,
           CTYPE1 = ctype1, CTYPE2 = ctype2,
           CUNIT1 = "deg",  CUNIT2 = "deg",
           CRVAL1 = RA,     CRVAL2 = Dec,
           CRPIX1 = (dim + 1) / 2, CRPIX2 = (dim + 1) / 2,
           CD1_1 = -ang * cos(t), CD1_2 =  ang * sin(t),
           CD2_1 =  ang * sin(t), CD2_2 =  ang * cos(t),
           NAXIS = 2L, NAXIS1 = as.integer(dim), NAXIS2 = as.integer(dim),
           RADESYS = "ICRS", EQUINOX = 2000),
      extra)
    class(kv) = "Rfits_keylist"
    kv
  }
  hdr <- function(kv) Rfits_header_to_raw(Rfits_keyvalues_to_header(kv))
  sipextra <- function(a20, a02, a11, a30 = 0) list(
    A_2_0 = a20, A_0_2 = a02, A_1_1 = a11, A_3_0 = a30,
    B_2_0 = a20 / 2, B_0_2 = a02 / 2, B_1_1 = -a11 / 2,
    APORDER = 3L, BPORDER = 3L)

  N <- 400L
  small_tan <- mk(RA = 180, Dec = 10, pixscale = 1, dim = N)

  cases <- list(
    "small TAN"          = list(hdr(small_tan), hdr(small_tan)),
    "small TAN shifted"  = list(hdr(small_tan),
                                hdr(mk(180.05, 10.02, 1, N))),
    "2x resample"        = list(hdr(small_tan), hdr(mk(180, 10, 0.5, N))),
    "rotated 45"         = list(hdr(mk(180, 10, 1, N, rot = 45)), hdr(small_tan)),
    "wide TAN 4.4deg"    = list(hdr(mk(180, 10, 20, N)), hdr(mk(180, 10, 20, N))),
    "wide TAN 10deg"     = list(hdr(mk(45, -30, 45, N)), hdr(mk(45, -30, 45, N))),
    "Dec 89.95"          = list(hdr(mk(180, 89.95, 20, N)), hdr(mk(180, 89.95, 20, N))),
    "Dec 89.999"         = list(hdr(mk(30, 89.999, 60, N)), hdr(mk(30, 89.999, 60, N))),
    "TPV bundled coeffs" = list(
      hdr(mk(180, 10, 20, N, ctype1 = "RA---ZPN", ctype2 = "DEC--ZPN",
             extra = list(PV2_1 = 1, PV2_2 = 0, PV2_3 = 44, PV2_5 = -10300))),
      hdr(mk(180, 10, 20, N, ctype1 = "RA---ZPN", ctype2 = "DEC--ZPN",
             extra = list(PV2_1 = 1, PV2_2 = 0, PV2_3 = 44, PV2_5 = -10300)))),
    "SIP cubic"          = list(
      hdr(mk(180, 10, 1, N, ctype1 = "RA---TAN-SIP", ctype2 = "DEC--TAN-SIP",
             extra = sipextra(1e-7, 2e-7, -1e-7, 1e-11))),
      hdr(mk(180, 10, 1, N, ctype1 = "RA---TAN-SIP", ctype2 = "DEC--TAN-SIP",
             extra = sipextra(1e-7, 2e-7, -1e-7, 1e-11)))),
    "SIP strong wide"    = list(
      hdr(mk(180, 10, 20, N, ctype1 = "RA---TAN-SIP", ctype2 = "DEC--TAN-SIP",
             extra = sipextra(5e-6, 8e-6, -6e-6, 1e-9))),
      hdr(mk(180, 10, 20, N, ctype1 = "RA---TAN-SIP", ctype2 = "DEC--TAN-SIP",
             extra = sipextra(5e-6, 8e-6, -6e-6, 1e-9)))),
    "wide+rot+SIP"       = list(
      hdr(mk(200, -45, 30, N, rot = 30, ctype1 = "RA---TAN-SIP",
             ctype2 = "DEC--TAN-SIP", extra = sipextra(2e-6, 3e-6, -2e-6, 5e-10))),
      hdr(mk(200, -45, 30, N, rot = 30, ctype1 = "RA---TAN-SIP",
             ctype2 = "DEC--TAN-SIP", extra = sipextra(2e-6, 3e-6, -2e-6, 5e-10)))),
    # Mismatched projection between in and out headers. These are expected to
    # refuse: the residual curvature across a cell is too large for bilinear.
    "TAN -> CAR"         = list(hdr(small_tan),
                                hdr(mk(180, 10, 20, N, ctype1 = "WC1-CAR",
                                       ctype2 = "WC2-CAR"))),
    "TAN -> SIN"         = list(hdr(small_tan),
                                hdr(mk(180, 10, 20, N, ctype1 = "WC1-SIN",
                                       ctype2 = "WC2-SIN")))
  )

  tol <- 1e-5
  accepted <- 0L
  for (lbl in names(cases)) {
    hin <- cases[[lbl]][[1]]; hout <- cases[[lbl]][[2]]
    wf <- function(x, y, cores = 1, ...)
      out2in(x, y, header_in = hin, header_out = hout, cores = cores)

    r <- suppressMessages(coarse(wf, c(N, N), tol = tol, step0 = 64L, cores = 1))

    if (is.null(r)) {
      # Refused: nothing to check, the caller falls back to the exact field.
      next
    }
    accepted <- accepted + 1L

    # The claimed error must be both true and inside the requested tolerance.
    pg <- expand.grid(seq_len(N), seq_len(N))
    ref <- suppressWarnings(wf(pg[, 1], pg[, 2], cores = 1))
    got <- cbind(as.vector(r$warpfield[, , 1]), as.vector(r$warpfield[, , 2]))
    ok <- is.finite(ref) & is.finite(got)
    true_err <- max(abs(got[ok] - ref[ok]))

    # expect_lt/expect_gte route `...` to compare(), not to an info argument, so
    # use expect_true to keep the failing geometry name in the report.
    expect_true(true_err < 1e-3, info = lbl)   # far better than a pixel
    expect_true(r$maxerr < tol, info = lbl)    # honoured its own contract
    # The midpoint estimate must not be wildly optimistic. A builder that
    # under-reported its error could pass tol while being wrong.
    expect_true(true_err < max(r$maxerr * 100, 1e-6), info = lbl)
  }

  # At least most geometries should be usable; if nearly everything refuses,
  # the coarse path is not earning its keep.
  expect_gte(accepted, 6L)
})

test_that("coarse warpgrid refuses when the transform is not finite", {
  suppressMessages(library(Rfits))
  ns <- asNamespace("ProPane")
  coarse <- get(".warpfield_coarse", envir = ns)
  out2in <- get(".warpfunc_out2in", envir = ns)

  src <- suppressMessages(Rfits_read_image(
    file.path(system.file("extdata/stack", package = "ProPane"), "image_1.fits")))
  ref <- suppressMessages(Rfits_read_image(
    file.path(system.file("extdata/stack", package = "ProPane"), "image_2.fits")))
  hi <- Rfits_keyvalues_to_raw(src$keyvalues)
  ho <- Rfits_header_to_raw(Rfits_keyvalues_to_header(ref$keyvalues))
  base <- function(x, y, cores = 1, ...)
    out2in(x, y, header_in = hi, header_out = ho, cores = cores)

  # A broad non-finite band stands in for a projection limb or a distortion
  # singularity: interpolation across it would invent finite warp values that
  # the exact path would instead have repaired with propanePatchPix.
  limb <- function(x, y, cores = 1, ...) {
    v <- base(x, y, cores = cores)
    hit <- y >= 150 & y <= 250
    v[hit, 1] <- Inf
    v
  }
  expect_null(suppressMessages(coarse(limb, c(400, 400), tol = 1e-5,
                                      step0 = 64L, cores = 1)))

  # Unreachable tolerance must refuse rather than silently under-deliver.
  expect_null(suppressMessages(
    coarse(base, c(400, 400), tol = 1e-30, step0 = 64L, max_refine = 2L,
           cores = 1)))
})
