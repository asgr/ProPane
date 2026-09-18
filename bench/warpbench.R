# propaneWarp benchmark harness.
#
# Reports min and median of n reps for each variant. Min is the headline number
# (least polluted by other work on the machine); median is shown too because
# this workload is allocation-heavy and the spread is real, not noise to ignore.
#
# NOTE for future editors: the timing helper must take a *function*. An earlier
# version took an expression and called force() inside the rep loop; promises
# memoise, so the work happened once and the result was divided by n, which
# reported absurdly fast stages and inverted the conclusions of the profiling.
#
#   Rscript bench/warpbench.R

local({
  suppressMessages(pkgload::load_all(".", quiet = TRUE, compile = TRUE,
                                     export_all = FALSE))
  suppressMessages({library(Rfits); library(Rwcs); library(imager)})

  ns <- asNamespace("ProPane")
  reps <- as.integer(Sys.getenv("BENCH_REPS", unset = "7"))

  d <- system.file("extdata/stack", package = "ProPane")
  frames <- sort(list.files(d, pattern = "^image_[0-9]+\\.fits$", full.names = TRUE))
  src <- suppressMessages(Rfits_read_image(frames[1]))
  ref <- suppressMessages(Rfits_read_image(frames[2]))
  kv <- ref$keyvalues
  dm <- dim(ref)[1:2]

  time1 <- function(f) {
    gc(FALSE)
    system.time(f())["elapsed"]
  }

  bench <- function(lbl, f) {
    invisible(f())  # warm
    t <- replicate(reps, time1(f))
    cat(sprintf("%-34s min %6.3fs  med %6.3fs  n=%d\n", lbl, min(t), median(t), reps))
    invisible(min(t))
  }

  P <- function(...) suppressMessages(propaneWarp(src, keyvalues_out = kv,
                                                   dim_out = dm, ...))

  cat(sprintf("input %s  output %dx%d   %s reps\n",
              paste(dim(src)[1:2], collapse = "x"), dm[1], dm[2], reps))
  cat("load avg:", paste(trimws(system("sysctl -n vm.loadavg", intern = TRUE)),
                         collapse = " "),
      " cores:", parallel::detectCores(), "\n\n")

  cat("=== single-band propaneWarp ===\n")
  base <- bench("warpgrid='exact' (default)", function() P())
  auto <- bench("warpgrid='auto'", function() P(warpgrid = "auto"))
  g64  <- bench("warpgrid=64", function() P(warpgrid = 64))
  nf   <- bench("warpgrid=64, dofinenorm=FALSE", function() P(warpgrid = 64, dofinenorm = FALSE))
  nc_  <- bench("warpgrid=64, dotightcrop=FALSE", function() P(warpgrid = 64, dotightcrop = FALSE))

  cat(sprintf("\n  auto    vs exact: %.2fx\n  grid64  vs exact: %.2fx\n",
              base / auto, base / g64))
  cat(sprintf("  finenorm costs: %.3fs  tightcrop costs: %.3fs\n",
              g64 - nf, g64 - nc_))

  cat("\n=== multi-band propaneWarpProPane ===\n")
  mkpp <- function() {
    x <- src
    structure(list(image = x, weight = x, inVar = x, exp = x, cold = x),
              class = "ProPane")
  }
  off <- bench("5 bands, share=FALSE",
               function() suppressMessages(propaneWarpProPane(
                 mkpp(), keyvalues_out = kv, dim_out = dm, warpfield_share = FALSE)))
  on  <- bench("5 bands, share=TRUE",
               function() suppressMessages(propaneWarpProPane(
                 mkpp(), keyvalues_out = kv, dim_out = dm, warpfield_share = TRUE)))
  onc <- bench("5 bands, share=TRUE, grid=64",
               function() suppressMessages(propaneWarpProPane(
                 mkpp(), keyvalues_out = kv, dim_out = dm, warpfield_share = TRUE,
                 warpgrid = 64)))
  cat(sprintf("\n  sharing:    %.2fx\n  + coarse:   %.2fx  (%.2fs -> %.2fs)\n",
              off / on, off / onc, off, onc))

  cat("\n=== 8-frame stack, 5 bands each ===\n")
  allf <- lapply(frames, function(f) suppressMessages(Rfits_read_image(f)))
  stack_bench <- function(share, warpgrid) {
    for (im in allf) {
      pp <- structure(list(image = im, weight = im, inVar = im, exp = im, cold = im),
                      class = "ProPane")
      suppressMessages(propaneWarpProPane(pp, keyvalues_out = kv, dim_out = dm,
                                          warpfield_share = share, warpgrid = warpgrid))
    }
    invisible(NULL)
  }
  s_off <- bench("8 frames, share=FALSE, exact", function() stack_bench(FALSE, "exact"))
  s_on  <- bench("8 frames, share=TRUE,  exact", function() stack_bench(TRUE, "exact"))
  s_c   <- bench("8 frames, share=TRUE,  grid=64", function() stack_bench(TRUE, 64))
  cat(sprintf("\n  total: %.2fx (sharing) -> %.2fx (sharing+coarse)  [%.1fs -> %.1fs]\n",
              s_off / s_on, s_off / s_c, s_off, s_c))
})
