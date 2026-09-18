# Generate propaneWarp golden fixtures from the CURRENT installed ProPane.
#
# Must be run BEFORE editing propaneWarp. Writes:
#   - tests/testthat/fixtures/warp-golden.rds : digests + samples (committed, small)
#   - /tmp/ProPane-warp-golden/<case>.rds     : full matrices (local only, for deep diffs)
#
# The committed digest is what makes the bit-identity tests exact; the local
# full data lets us measure the magnitude of a deliberate difference.

local({
  pkg <- "ProPane"
  suppressMessages({
    library(package = pkg, character.only = TRUE)
    library(Rfits); library(Rwcs); library(imager)
  })

  ns <- asNamespace(pkg)
  propaneWarp <- get("propaneWarp", envir = ns)

  # Hash the raw double bytes rather than a serialised object, so the digest is
  # independent of R's serialization format, and use tools::md5sum (base R)
  # rather than digest: CI installs hard dependencies only, so anything in
  # Suggests may be absent when the tests run.
  digest_matrix <- function(m) {
    v <- as.vector(m, mode = "double")
    f <- tempfile()
    con <- file(f, "wb")
    # Native byte order, fixed width: keeps NaN bit patterns intact.
    writeBin(v, con, size = 8)
    close(con)
    on.exit(unlink(f))
    unname(tools::md5sum(f))
  }

  stackdir <- system.file("extdata/stack", package = pkg)
  frames <- sort(list.files(stackdir, pattern = "^image_[0-9]+\\.fits$", full.names = TRUE))
  stopifnot(length(frames) >= 2)

  im_src <- suppressMessages(Rfits_read_image(frames[1]))
  im_ref <- suppressMessages(Rfits_read_image(frames[2]))
  dim_ref <- dim(im_ref)[1:2]

  # Deliberately-off target WCS offsets, used for the degenerate-crop probes
  # below (see `probe`).

  # Cases whose 1.10.1 output is NOT reproducible run-to-run.
  # direction="forward" uses CImg's forward scatter, which accumulates into the
  # output without atomics: repeat runs at default threads give different
  # results, and the answer also depends on OMP_NUM_THREADS. Backward warp is a
  # pure gather and is bit-stable, so only the explicit-forward case is here
  # ("auto" resolves to backward whenever pixscale_in >= pixscale_out).
  NONDETERMINISTIC <- "forward"
  ND_REL_TOL <- 1e-3

  cases <- list(
    default        = list(),
    forward        = list(direction = "forward"),
    nearest        = list(interpolation = "nearest"),
    linear         = list(interpolation = "linear"),
    neumann        = list(boundary = "neumann"),
    periodic       = list(boundary = "periodic"),
    nofinenorm     = list(dofinenorm = FALSE),
    notightcrop    = list(dotightcrop = FALSE),
    keepcrop       = list(keepcrop = TRUE),
    extratight     = list(keepcrop = TRUE, extratight = TRUE),
    noscale        = list(doscale = FALSE),
    magzero        = list(magzero_out = 25, magzero_in = 30),
    warpfield_in   = list(prebuild_field = TRUE)
  )

  outdir_full <- "/tmp/ProPane-warp-golden"
  dir.create(outdir_full, showWarnings = FALSE, recursive = TRUE)

  fixt <- list()

  for (nm in names(cases)) {
    arg <- cases[[nm]]
    prebuild <- isTRUE(arg$prebuild_field)
    arg$prebuild_field <- NULL
    kvo <- im_ref$keyvalues

    t0 <- proc.time()[3]
    if (prebuild) {
      arg$warpfield <- suppressMessages(
        propaneWarp(im_src, keyvalues_out = kvo, dim_out = dim_ref,
                    warpfield_return = TRUE))$warpfield
    }

    r <- tryCatch(
      suppressMessages(do.call(propaneWarp,
        c(list(image_in = im_src, keyvalues_out = kvo, dim_out = dim_ref), arg))),
      error = function(e) e
    )
    dt <- proc.time()[3] - t0

    if (inherits(r, "try-error") || inherits(r, "error")) {
      fixt[[nm]] <- list(status = "error", msg = conditionMessage(r))
      cat(sprintf("%-14s ERROR  %s\n", nm, conditionMessage(r)))
      next
    }

    mat <- r$imDat
    saveRDS(mat, file.path(outdir_full, paste0(nm, ".rds")))

    # Re-run to find out, empirically, whether this case is reproducible.
    d0 <- digest_matrix(mat)
    digs <- replicate(2, digest_matrix(
      suppressMessages(do.call(propaneWarp,
        c(list(image_in = im_src, keyvalues_out = kvo, dim_out = dim_ref), arg)))$imDat))
    repro <- all(digs == d0)
    if (!repro) {
      cat(sprintf("%-14s NON-DETERMINISTIC in 1.10.1 (%d distinct digests over 3 runs)\n",
                  nm, length(unique(c(digs, d0)))))
    }

    set.seed(20260916L)
    idx <- sort(sample.int(length(mat), min(4096L, length(mat))))
    fixt[[nm]] <- list(
      status = if (repro) "ok" else "nondeterministic",
      reproducible = repro,
      rel_tol = if (repro) 0 else ND_REL_TOL,
      digest = d0,
      dim = dim(mat),
      sum = sum(mat, na.rm = TRUE),
      nNA = sum(is.na(mat)),
      sample_idx = idx,
      sample_val = mat[idx],
      kv = r$keyvalues[intersect(
        c("NAXIS1","NAXIS2","CRPIX1","CRPIX2","XCUTLO","XCUTHI","YCUTLO","YCUTHI","MAGZERO"),
        names(r$keyvalues))]
    )
    cat(sprintf("%-14s %s  %-12s %.2fs  sum=%.10g nNA=%d\n", nm,
                if (repro) "ok  " else "NONDET",
                paste(dim(mat), collapse = "x"), dt, fixt[[nm]]$sum, fixt[[nm]]$nNA))
  }

  # Degenerate crop: in 1.10.1 an inverted tight-crop range is passed straight
  # to the box crop, which silently *expands* the internal working grid. The
  # returned image still looks normal-sized, so we record the warpfield dims
  # (which expose the internal grid) alongside the output digest.
  #
  # These entries are a HISTORICAL RECORD of the 1.10.1 blowup, not a target:
  # propaneWarp no longer warps a non-overlapping frame at all (it returns the
  # requested frame filled with `blank`), so the test asserts on out_dim and
  # the header keys -- which should still match -- and treats digest, field_dim
  # and field_px as evidence of what used to happen.
  probe <- function(sh) {
    kv <- im_ref$keyvalues
    kv$CRPIX1 <- kv$CRPIX1 + sh
    t0 <- proc.time()[3]
    m <- suppressMessages(propaneWarp(im_src, keyvalues_out = kv, dim_out = dim_ref,
                                      warpfield_return = TRUE))
    list(shift = sh,
         status = "ok",
         field_dim = dim(m$warpfield)[1:2],
         field_px = prod(dim(m$warpfield)[1:2]),
         out_dim = dim(m$imDat),
         digest = digest_matrix(m$imDat),
         nNA = sum(is.na(m$imDat)),
         allNA = all(is.na(m$imDat)),
         sum = sum(m$imDat, na.rm = TRUE),
         kv = m$keyvalues[c("NAXIS1","NAXIS2","CRPIX1","CRPIX2",
                            "XCUTLO","XCUTHI","YCUTLO","YCUTHI")],
         wall = proc.time()[3] - t0)
  }

  for (sh in c(2000L, -10000L)) {
    p <- probe(sh)
    key <- paste0("degenerate_", ifelse(sh > 0, "pos", "neg"))
    fixt[[key]] <- p
    cat(sprintf("%-15s field %-9s (%s px)  out %s  allNA=%s  %.2fs\n", key,
                paste(p$field_dim, collapse = "x"), format(p$field_px, big.mark = ","),
                paste(p$out_dim, collapse = "x"), p$allNA, p$wall))
  }
  # Also record the well-behaved reference warpfield size for the same geometry.
  fixt[["degenerate_reference_field_px"]] <-
    prod(dim(suppressMessages(propaneWarp(im_src, keyvalues_out = im_ref$keyvalues,
                                          dim_out = dim_ref,
                                          warpfield_return = TRUE))$warpfield)[1:2])
  cat(sprintf("reference field px: %s\n",
              format(fixt[["degenerate_reference_field_px"]], big.mark = ",")))

  retdir <- file.path("tests", "testthat", "fixtures")
  dir.create(retdir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(list(
    created = Sys.time(),
    propaneversion = as.character(packageVersion(pkg)),
    src_frame = basename(frames[1]), ref_frame = basename(frames[2]),
    dim_ref = dim_ref, cases = fixt
  ), file.path(retdir, "warp-golden.rds"), version = 2)
  cat("\nwrote", file.path(retdir, "warp-golden.rds"), "\n")
})
