# Accuracy gate for the coarse warp-field path (.warpfield_coarse).
#
# The bundled stack frames are a small 10x10 arcmin near-affine field, which
# almost any approximation passes. This exercises WCS types that stress the
# bilinear reconstruction: large spherical curvature, strong distortion, high
# declination, steep rotation, and mismatched projections.
#
# Headers are synthesised in memory via propaneGenWCS, so nothing needs to be
# shipped in inst/extdata. Run with:
#   Rscript bench/stress-warpgrid.R

local({
  pkg <- "ProPane"
  # Load from source so this script tracks the working tree, not the installed
  # build (the coarse path is internal and only exists in-tree for now).
  suppressMessages(pkgload::load_all(".", quiet = TRUE, compile = TRUE,
                                     export_all = FALSE))
  suppressMessages({
    library(Rfits); library(Rwcs); library(imager)
  })
  ns <- asNamespace(pkg)
  warpfield_coarse <- get(".warpfield_coarse", envir = ns)
  warpfunc_out2in  <- get(".warpfunc_out2in", envir = ns)

  # Compare the coarse reconstruction against the full-resolution field, which
  # is the thing the caller would otherwise have paid for.
  check <- function(lbl, header_in, header_out, dim_xy, tol = 1e-5, step0 = 64L) {
    nx <- dim_xy[1]; ny <- dim_xy[2]
    wf <- function(x, y, cores = 1, ...)
      warpfunc_out2in(x, y, header_in = header_in, header_out = header_out,
                      cores = cores)

    r <- tryCatch({
      t0 <- proc.time()[3]
      pg <- expand.grid(seq_len(nx), seq_len(ny))
      ref <- suppressWarnings(wf(pg[, 1], pg[, 2], cores = 1))
      t_exact <- proc.time()[3] - t0

      t1 <- proc.time()[3]
      got <- suppressMessages(tryCatch(
        warpfield_coarse(wf, c(nx, ny), tol = tol, step0 = step0, cores = 1),
        error = function(e) structure(conditionMessage(e), class = "err")))
      t_coarse <- proc.time()[3] - t1
      list(ref = ref, got = got, t_exact = t_exact, t_coarse = t_coarse)
    }, error = function(e) {
      cat(sprintf("%-30s  HEADER/WCS ERROR: %s\n", lbl,
                  gsub("\n", " ", conditionMessage(e))))
      NULL
    })

    if (is.null(r)) return(invisible(NULL))
    if (inherits(r$got, "err")) {
      cat(sprintf("%-30s  ERROR: %s\n", lbl, gsub("\n", " ", r$got)))
      return(invisible(NULL))
    }
    if (is.null(r$got)) {
      cat(sprintf("%-30s  REFUSED -> falls back to exact (%.2fs)\n", lbl, r$t_exact))
      return(invisible(NULL))
    }

    got_x <- as.vector(r$got$warpfield[, , 1])
    got_y <- as.vector(r$got$warpfield[, , 2])
    ok <- is.finite(r$ref[, 1]) & is.finite(r$ref[, 2]) &
          is.finite(got_x)      & is.finite(got_y)
    err <- max(abs(got_x[ok] - r$ref[ok, 1]), abs(got_y[ok] - r$ref[ok, 2]))

    cat(sprintf("%-30s step %3d  maxerr %9.2e px  %-5s  exact %.2fs -> coarse %.3fs (%.1fx)\n",
                lbl, r$got$step, err,
                if (err < 0.01) "PASS" else if (err < 0.1) "WARN" else "FAIL",
                r$t_exact, r$t_coarse, r$t_exact / r$t_coarse))
    invisible(list(label = lbl, step = r$got$step, err = err,
                   claimed = r$got$maxerr))
  }

  # Synthesise a WCS keyvalue list from scratch. propaneGenWCS() insists on
  # scanning real files, which is the wrong dependency for a synthetic test.
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

  N <- 800L

  cat("=== baseline ===\n")
  small_tan <- mk(RA = 180, Dec = 10, pixscale = 1, dim = N)
  check("small TAN (0.22 deg)", hdr(small_tan), hdr(small_tan), c(N, N))

  shift <- mk(RA = 180.05, Dec = 10.02, pixscale = 1, dim = N)
  check("small TAN, shifted", hdr(small_tan), hdr(shift), c(N, N))

  scale <- mk(RA = 180, Dec = 10, pixscale = 0.5, dim = N)
  check("TAN 2x resample", hdr(small_tan), hdr(scale), c(N, N))

  cat("\n=== curvature / size stress ===\n")
  wide <- mk(RA = 180, Dec = 10, pixscale = 20, dim = N)   # ~4.4 deg
  check("wide TAN (4.4 deg)", hdr(wide), hdr(wide), c(N, N))

  wide2 <- mk(RA = 45, Dec = -30, pixscale = 45, dim = N)  # ~10 deg
  check("very wide TAN (10 deg)", hdr(wide2), hdr(wide2), c(N, N))

  rot <- mk(RA = 180, Dec = 10, pixscale = 1, dim = N, rot = 45)
  check("TAN rotated 45 deg", hdr(rot), hdr(small_tan), c(N, N))

  cat("\n=== projection stress ===\n")
  # A single header cannot mix projections (wcslib rejects CTYPE1=TAN with
  # CTYPE2=CAR), so the mismatch is between the in and out WCS.
  car <- mk(RA = 180, Dec = 10, pixscale = 20, dim = N,
            ctype1 = "WC1-CAR", ctype2 = "WC2-CAR")
  check("TAN -> CAR (mismatched)", hdr(small_tan), hdr(car), c(N, N))

  sin_ <- mk(RA = 180, Dec = 10, pixscale = 20, dim = N,
             ctype1 = "WC1-SIN", ctype2 = "WC2-SIN")
  check("TAN -> SIN (mismatched)", hdr(small_tan), hdr(sin_), c(N, N))

  gp <- mk(RA = 180, Dec = 10, pixscale = 20, dim = N,
           ctype1 = "RA---GLQ", ctype2 = "DEC--GLQ")
  check("TAN -> GLQ (mismatched)", hdr(small_tan), hdr(gp), c(N, N))

  cat("\n=== near pole ===\n")
  pole <- mk(RA = 180, Dec = 89.95, pixscale = 20, dim = N)
  check("TAN at Dec 89.95", hdr(pole), hdr(pole), c(N, N))

  pole2 <- mk(RA = 30, Dec = 89.999, pixscale = 60, dim = N)
  check("TAN at Dec 89.999", hdr(pole2), hdr(pole2), c(N, N))

  cat("\n=== distortion stress ===\n")
  tpv <- mk(RA = 180, Dec = 10, pixscale = 20, dim = N,
            ctype1 = "RA---ZPN", ctype2 = "DEC--ZPN",
            extra = list(PV2_1 = 1, PV2_2 = 0, PV2_3 = 44, PV2_5 = -10300))
  check("TPV (bundled coefficients)", hdr(tpv), hdr(tpv), c(N, N))

  tpv2 <- mk(RA = 45, Dec = -30, pixscale = 45, dim = N,
             ctype1 = "RA---TPV", ctype2 = "DEC--TPV",
             extra = list(PV1_2 = 30, PV1_3 = -500, PV2_2 = 25, PV2_3 = -400))
  check("strong TPV on wide field", hdr(tpv2), hdr(tpv2), c(N, N))

  sipextra <- function(a20, a02, a11, a30 = 0, aporder = 3) list(
    A_2_0 = a20, A_0_2 = a02, A_1_1 = a11, A_3_0 = a30,
    B_2_0 = a20 / 2, B_0_2 = a02 / 2, B_1_1 = -a11 / 2,
    APORDER = as.integer(aporder), BPORDER = as.integer(aporder))

  sip <- mk(RA = 180, Dec = 10, pixscale = 1, dim = N,
            ctype1 = "RA---TAN-SIP", ctype2 = "DEC--TAN-SIP",
            extra = sipextra(1e-7, 2e-7, -1e-7, 1e-11))
  check("SIP cubic distortion", hdr(sip), hdr(sip), c(N, N))

  sip2 <- mk(RA = 180, Dec = 10, pixscale = 20, dim = N,
             ctype1 = "RA---TAN-SIP", ctype2 = "DEC--TAN-SIP",
             extra = sipextra(5e-6, 8e-6, -6e-6, 1e-9))
  check("strong SIP on wide field", hdr(sip2), hdr(sip2), c(N, N))

  cat("\n=== combined: wide + rotated + distorted ===\n")
  hard <- mk(RA = 200, Dec = -45, pixscale = 30, dim = N, rot = 30,
             ctype1 = "RA---TAN-SIP", ctype2 = "DEC--TAN-SIP",
             extra = sipextra(2e-6, 3e-6, -2e-6, 5e-10))
  check("wide+rotated+SIP", hdr(hard), hdr(hard), c(N, N))
})
