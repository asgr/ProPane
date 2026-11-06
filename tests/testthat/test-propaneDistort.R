# Pixel-scale / pixel-area maps.
#
# Two properties matter and both used to fail on a field crossing RA = 0:
#
# 1. The answer must not depend on where the field sits in longitude. The old
#    code multiplied each corner's *absolute* RA by cos(dec), so shifting the
#    same field by 180 deg changed the reported area by a factor of 2.37 (on the
#    WWS mosaic that prompted this).
# 2. RA differences must be wrapped across the 0/360 seam. The old code took
#    BR[,1] - BL[,1] literally, so pixels straddling the seam reported areas of
#    up to 3.3e8 asec2 where the true maximum is 9.0e4.
#
# Rfits::pixarea() / pixscale() are the independent reference: they already wrap
# the seam and use the parallelogram (cross-product) area. They are hard
# dependencies, so no guarded skip is needed.
#
# Headers are synthesised rather than read from inst/extdata, so this file costs
# a fraction of a second and stays free of foreach/Rfits_make_list (see the
# header of test-propaneStackWarpInVar.R for why that matters in this session).

.distort_kv = function(RA, Dec = -37.5, pixscale = 300, dim = c(120L, 60L)) {
  suppressMessages(library(Rfits))
  ang = pixscale/3600
  kv = list(CTYPE1 = 'RA---TAN', CTYPE2 = 'DEC--TAN',
            CUNIT1 = 'deg', CUNIT2 = 'deg',
            CRVAL1 = RA, CRVAL2 = Dec,
            CRPIX1 = (dim[1] + 1)/2, CRPIX2 = (dim[2] + 1)/2,
            CD1_1 = -ang, CD1_2 = 0, CD2_1 = 0, CD2_2 = ang,
            NAXIS = 2L, NAXIS1 = dim[1], NAXIS2 = dim[2],
            RADESYS = 'ICRS', EQUINOX = 2000)
  class(kv) = 'Rfits_keylist'
  kv
}

# Sampled pixels, corners included, so the extreme distortion is covered.
.probe_pix = function(dim) {
  expand.grid(x = c(1L, round(dim[1]/3), round(dim[1]/2), dim[1]),
              y = c(1L, round(dim[2]/2), dim[2]))
}

test_that("pixarea and pixscale are invariant to a longitude shift", {
  skip_on_cran()
  kv_ref = .distort_kv(RA = 160)    # 20 deg west of this field's 180 deg meridian
  kv_shift = .distort_kv(RA = -20)  # the same field, 180 deg away
  # RA = -20 and 340 describe the same sky position; keep the shiftee inside
  # [0, 360) so the header is legal.
  kv_shift$CRVAL1 = 340

  area_ref = propaneDistortPixarea(kv_ref)$imDat
  area_shift = propaneDistortPixarea(kv_shift)$imDat
  scale_ref = propaneDistortPixscale(kv_ref)$imDat
  scale_shift = propaneDistortPixscale(kv_shift)$imDat

  expect_equal(max(abs(area_shift/area_ref - 1)), 0, tolerance = 1e-9)
  expect_equal(max(abs(scale_shift/scale_ref - 1)), 0, tolerance = 1e-9)
})

test_that("the RA = 0 seam does not inflate the maps", {
  skip_on_cran()
  # 10 deg wide centred on RA 3, so the seam runs through the middle of the field
  kv_seam = .distort_kv(RA = 3)
  kv_away = .distort_kv(RA = 183)

  area_seam = propaneDistortPixarea(kv_seam)$imDat
  area_away = propaneDistortPixarea(kv_away)$imDat

  # A seam crossing must be a no-op: peak-to-peak range and every value agree
  # with the identical field placed away from the seam.
  expect_equal(max(abs(area_seam/area_away - 1)), 0, tolerance = 1e-9)
  # Before the fix the seam column was ~3600x the neighbouring column. The real
  # field varies by well under an order of magnitude, so a loose bound catches a
  # regression without pinning the exact geometry.
  expect_lt(max(area_seam)/min(area_seam), 10)
})

test_that("pixarea and pixscale agree with Rfits at sampled pixels", {
  skip_on_cran()
  for (RA in c(3, 183)) {          # seam-crossing and not
    kv = .distort_kv(RA = RA)
    img = list(keyvalues = kv)
    class(img) = 'Rfits_image'

    area = propaneDistortPixarea(kv)$imDat
    scale = propaneDistortPixscale(kv)$imDat

    pts = .probe_pix(dim(kv))
    # Rfits loc is in the 0-based pixel-centre convention: the centre of FITS
    # pixel i is loc = i - 0.5.
    ref_area = mapply(function(x, y) suppressMessages(Rfits::pixarea(img, loc = c(x, y))),
                      pts$x - 0.5, pts$y - 0.5)
    ref_scale = mapply(function(x, y) suppressMessages(Rfits::pixscale(img, loc = c(x, y))),
                       pts$x - 0.5, pts$y - 0.5)

    ours_area = mapply(function(x, y) area[x, y], pts$x, pts$y)
    ours_scale = mapply(function(x, y) scale[x, y], pts$x, pts$y)

    expect_equal(max(abs(ours_area/ref_area - 1)), 0, tolerance = 1e-9)
    expect_equal(max(abs(ours_scale/ref_scale - 1)), 0, tolerance = 1e-9)
  }
})

test_that("units are consistent and the output matches the input grid", {
  skip_on_cran()
  kv = .distort_kv(RA = 3)

  asec2 = propaneDistortPixarea(kv, unit = 'asec2')$imDat
  deg2 = propaneDistortPixarea(kv, unit = 'deg2')$imDat
  amin2 = propaneDistortPixarea(kv, unit = 'amin2')$imDat
  str = propaneDistortPixarea(kv, unit = 'str')$imDat

  expect_s3_class(propaneDistortPixarea(kv), 'Rfits_image')
  expect_identical(dim(asec2), dim(kv))
  expect_equal(deg2*3600^2, asec2, tolerance = 1e-12)
  expect_equal(amin2*60^2, asec2, tolerance = 1e-12)
  expect_equal(str, deg2*(pi/180)^2, tolerance = 1e-12)

  asec = propaneDistortPixscale(kv, unit = 'asec')$imDat
  amin = propaneDistortPixscale(kv, unit = 'amin')$imDat
  expect_equal(amin*60, asec, tolerance = 1e-12)
  # A scale is the side of the equal-length square, so its square is the area
  # only for a square pixel; they must still be the same order here.
  expect_lt(max(abs(sqrt(asec2)/asec - 1)), 0.05)
})