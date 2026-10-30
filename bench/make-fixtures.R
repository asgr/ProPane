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

  # Cases that must ALWAYS be recorded as non-reproducible, whatever the local
  # re-run below happens to observe. Empty at present; see below for why.
  #
  # HISTORY, because this bit us hard and the reasoning still matters.
  #
  # Until imager 1.0.10, `direction = "forward"` was here. It uses CImg's
  # forward-absolute scatter, which blended into the output
  # (set_linear_atXY) with a read-modify-write that had no atomics. Whether that
  # was parallelised was a property of the CImg build, not of this package:
  # stock CImg wrapped the scatter in cimg_pragma_openmp(parallel for ...) so
  # contributions raced and were lost, while the asgr fork patch serialised it
  # and hid the problem entirely.
  #
  # The trap was that measuring reproducibility by re-running inside the same
  # process can only ever test the build it runs on. On a serial build the
  # re-run always agreed, the case got recorded as exact and gated at
  # GOLDEN_TOL (1e-9), and it then failed on every stock-CImg machine -- which
  # is precisely what turned CI red. Measured cross-build spread then: flux
  # agreed to 2.9e-5 and ~1.6% of sampled pixels differed (max 3.4e-2
  # relative), a racy lost update rather than rounding.
  #
  # imager 1.0.10 (CImg 4.0.7) removes the reason for the exemption: the scatter
  # now accumulates with is_added into a separate `accu` count and normalises
  # once at the end, in a fixed serial order. That is deterministic AND correct
  # (the old blend degenerated to last-write-wins at integer destinations), so
  # forward is now measured like every other case and pinned at GOLDEN_TOL.
  #
  # Keep this mechanism rather than deleting it: if a future CImg parallelises
  # the accumulation, or reintroduces an order-dependent blend, the fix is to
  # list the case here again -- not to widen a tolerance until CI goes green.
  NONDETERMINISTIC <- character(0)
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
    warpfield_in   = list(prebuild_field = TRUE),
    # Contracting geometry (~4 sources per destination), the regime forward
    # warping exists for. The loop below swaps in contraction_kv/dim_contraction
    # for this name; the args are the same forward direction as `forward`.
    contraction_fwd = list(direction = "forward")
  )

  outdir_full <- "/tmp/ProPane-warp-golden"
  dir.create(outdir_full, showWarnings = FALSE, recursive = TRUE)

  fixt <- list()

  # ---------------------------------------------------------------------------
  # Contraction geometry, used for the `contraction_fwd` case below.
  #
  # Every other case warps onto the reference frame, whose pixel scale is within
  # 3e-5 of the input's (0.3372838 vs 0.3372744). That is a near-identity
  # resampling, so at most a couple of source pixels share a destination and the
  # forward scatter's behaviour is barely exercised. The forward-warp bug fixed
  # in imager 1.0.10 -- last-write-wins instead of averaging -- therefore shows
  # up there as a small pixel-level change with only a 2.6e-5 flux difference.
  #
  # This case doubles the output pixel scale instead, so roughly four source
  # pixels land in every destination. That is the regime forward warping exists
  # for (drizzling onto a coarser grid), and the regime where getting the
  # scatter wrong actually loses data. It is derived from the reference header
  # so it shares its sky position: CD is doubled, the pixel size halves in each
  # axis, and CRPIX is adjusted to keep the same sky point at the frame centre.
  # ---------------------------------------------------------------------------
  contraction_kv <- im_ref$keyvalues
  contraction_kv$CD1_1 <- im_ref$keyvalues$CD1_1 * 2
  contraction_kv$CD1_2 <- im_ref$keyvalues$CD1_2 * 2
  contraction_kv$CD2_1 <- im_ref$keyvalues$CD2_1 * 2
  contraction_kv$CD2_2 <- im_ref$keyvalues$CD2_2 * 2
  contraction_kv$CRPIX1 <- (im_ref$keyvalues$CRPIX1 + 0.5) / 2 - 0.5
  contraction_kv$CRPIX2 <- (im_ref$keyvalues$CRPIX2 + 0.5) / 2 - 0.5
  dim_contraction <- ceiling(dim_ref / 2)

  for (nm in names(cases)) {
    arg <- cases[[nm]]
    prebuild <- isTRUE(arg$prebuild_field)
    arg$prebuild_field <- NULL
    kvo <- im_ref$keyvalues
    dim_use <- dim_ref

    if (nm == "contraction_fwd") {
      kvo <- contraction_kv
      dim_use <- dim_contraction
    }

    t0 <- proc.time()[3]
    if (prebuild) {
      arg$warpfield <- suppressMessages(
        propaneWarp(im_src, keyvalues_out = kvo, dim_out = dim_use,
                    warpfield_return = TRUE))$warpfield
    }

    r <- tryCatch(
      suppressMessages(do.call(propaneWarp,
        c(list(image_in = im_src, keyvalues_out = kvo, dim_out = dim_use), arg))),
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

    # Find out empirically whether this case is reproducible. That measurement is
    # only trustworthy for a case whose determinism is a property of the
    # algorithm rather than of the build, which is now true of every case here
    # (see the NONDETERMINISTIC note above); anything listed there would have
    # its answer overridden below.
    d0 <- digest_matrix(mat)
    digs <- replicate(2, digest_matrix(
      suppressMessages(do.call(propaneWarp,
        c(list(image_in = im_src, keyvalues_out = kvo, dim_out = dim_use), arg)))$imDat))
    repro <- all(digs == d0)
    if (!repro) {
      cat(sprintf("%-14s NON-DETERMINISTIC on this build (%d distinct digests over 3 runs)\n",
                  nm, length(unique(c(digs, d0)))))
    }
    if (nm %in% NONDETERMINISTIC) {
      if (repro) {
        cat(sprintf(paste0(
          "%-14s build-dependent path, but this build agreed over 3 runs; recording\n",
          "               it as nondeterministic anyway.\n",
          "               Pinning it here would fail on a build where it is racy.\n"), nm))
      }
      repro <- FALSE
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
      # ZNAXIS1/2 are recorded alongside NAXIS1/2 because these frames are
      # tile-compressed: on a ZIMAGE frame NAXIS1/2 are the BINTABLE compressed
      # row width (bytes) and row count, while the image size lives in ZNAXIS1/2.
      # Recording both lets the test assert on whichever pair is authoritative
      # without special-casing, and makes the distinction visible in the fixture.
      kv = r$keyvalues[intersect(
        c("ZIMAGE", "NAXIS1","NAXIS2","ZNAXIS1","ZNAXIS2",
          "CRPIX1","CRPIX2","XCUTLO","XCUTHI","YCUTLO","YCUTHI","MAGZERO"),
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
         kv = m$keyvalues[intersect(
           c("ZIMAGE", "NAXIS1","NAXIS2","ZNAXIS1","ZNAXIS2",
             "CRPIX1","CRPIX2","XCUTLO","XCUTHI","YCUTLO","YCUTHI"),
           names(m$keyvalues))],
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
    dim_ref = dim_ref,
    # Stored rather than re-derived in the test, so the fixture and the test
    # cannot drift apart on what "the contraction geometry" means.
    contraction_kv = contraction_kv,
    dim_contraction = dim_contraction,
    cases = fixt
  ), file.path(retdir, "warp-golden.rds"), version = 2)
  cat("\nwrote", file.path(retdir, "warp-golden.rds"), "\n")
})
