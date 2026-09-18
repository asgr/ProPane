.warpfunc_in2out = function(x, y, header_in=NULL, WCSref_in=NULL, header_out=NULL, WCSref_out=NULL, cores=1) {
  radectemp = Rwcs_p2s(x, y, header = header_in, WCSref = WCSref_in, cores = cores)
  xy_out = Rwcs_s2p(radectemp[,1], radectemp[,2], header = header_out, WCSref = WCSref_out, cores = cores)
  return(xy_out)
}
.warpfunc_out2in = function(x, y, header_in=NULL, WCSref_in=NULL, header_out=NULL, WCSref_out=NULL, cores=1) {
  radectemp = Rwcs_p2s(x, y, header = header_out, WCSref = WCSref_in, cores = cores)
  xy_out = Rwcs_s2p(radectemp[,1], radectemp[,2], header = header_in, WCSref = WCSref_out, cores = cores)
  return(xy_out)
}

# ---------------------------------------------------------------------------
# Coarse-grid warp field construction.
#
# propaneWarp() normally evaluates the celestial transform at *every* output
# pixel, which is by far the dominant cost: measured on a 1816x1767 output the
# two Rwcs calls are 1.55s of a 2.16s call (~235 ns/pixel, inside wcslib).
#
# The out->in mapping is a smooth function of position, so it can be recovered
# to well below sub-pixel accuracy from a sparse lattice:
#
#   field(i, j) = (a0 + a1*i + a2*j)   least-squares affine, exact to first order
#               + r(i, j)              bilinear interpolation of the residual
#
# Subtracting the affine before interpolating is what makes this accurate: the
# residual is left with only the field curvature, so bilinear error is tiny even
# on a coarse lattice. Accuracy is verified a posteriori against the true WCS at
# cell midpoints, and we refine (or give up) rather than trust a heuristic.
#
# Returns a list(warpfield, step, maxerr), or NULL when the field cannot be
# built to tolerance -- callers must then fall back to the exact path.
.warpfield_coarse = function(warpfun, dim_xy, tol = 1e-5, step0 = 64, cores = 1,
                             max_refine = 5, ...){
  nx = dim_xy[1]; ny = dim_xy[2]
  if(nx < 4L || ny < 4L) return(NULL)

  fit_affine = function(gx, gy, gv){
    A = cbind(1, gx, gy)
    b = qr.solve(crossprod(A), crossprod(A, gv))
    list(b = as.numeric(b), r = gv - as.numeric(A %*% b))
  }

  # A-posteriori verification, in two parts.
  #
  # (1) Accuracy. Lattice nodes are reproduced exactly by construction, so error
  # only accumulates inside a cell. Every cell interior is probed at
  # (0.25, 0.5, 0.75) along each axis -- 9 points per cell -- and the
  # reconstruction is compared with the true transform there. Measured against
  # the full-resolution field this estimate is accurate to within ~1%, so no
  # fudge factor is needed.
  #
  # (2) Finiteness. A region where the true transform is not finite is a hard
  # failure: the exact path repairs those pixels spatially with propanePatchPix,
  # which a sparse lattice cannot reproduce, and interpolating across one would
  # invent finite warp values. Accuracy probes alone can miss a narrow band --
  # at step 64 they are 16px apart, so a 3px band slips through. So we sweep a
  # separate finiteness grid at a fixed spacing (LIMB_CHECK_STEP px), which caps
  # the width of any undetected contiguous non-finite region at that spacing.
  # Real projection limbs and distortion singularities are broad, but a small
  # independent check is cheap: at 8px this is nx*ny/64 points, ~50k, about
  # 0.03s, versus 1.5s for the full per-pixel transform. It is a bound on the
  # failure width, not a proof of absence.
  FRACS = c(0.25, 0.5, 0.75)
  LIMB_CHECK_STEP = 8L

  s = max(2L, min(as.integer(step0), max(2L, min(nx, ny) %/% 2L)))

  for(attempt in seq_len(max_refine + 1L)){
    cxs = seq.int(1L, nx, by = s)
    cys = seq.int(1L, ny, by = s)
    if(cxs[length(cxs)] != nx) cxs = c(cxs, nx)
    if(cys[length(cys)] != ny) cys = c(cys, ny)
    nc = length(cxs); nr = length(cys)
    if(nc < 3L || nr < 3L) return(NULL)

    gx = rep(cxs, times = nr)
    gy = rep(cys, each = nc)
    co = warpfun(gx, gy, cores = cores, ...)

    # Infinities are exactly what the exact path repairs spatially with
    # propanePatchPix, which cannot be reproduced from a sparse lattice.
    if(!all(is.finite(co))) return(NULL)

    # Dense finiteness sweep -- see part (2) of the verification note above.
    # Only needs to run once per call, and is independent of the lattice step.
    if(attempt == 1L){
      fxv = seq.int(1L, nx, by = LIMB_CHECK_STEP)
      fyv = seq.int(1L, ny, by = LIMB_CHECK_STEP)
      sweep = suppressWarnings(warpfun(rep(fxv, times = length(fyv)),
                                       rep(fyv, each = length(fxv)),
                                       cores = cores, ...))
      if(!all(is.finite(sweep))) return(NULL)
    }

    fx = fit_affine(gx, gy, co[, 1])
    fy = fit_affine(gx, gy, co[, 2])

    # cell-relative probe coordinates; entry [a, k] = cxs[a] + dx[a]*FRACS[k],
    # flattened so that a varies slowest and k fastest.
    xpr = as.vector(t(cxs[-nc] + outer(diff(cxs), FRACS, FUN = '*')))
    ypr = as.vector(t(cys[-nr] + outer(diff(cys), FRACS, FUN = '*')))
    px = rep(xpr, times = length(ypr))
    py = rep(ypr, each = length(xpr))

    po = suppressWarnings(warpfun(px, py, cores = cores, ...))
    if(!all(is.finite(po))) return(NULL)

    cellx = pmin(findInterval(px, cxs), nc - 1L)
    celly = pmin(findInterval(py, cys), nr - 1L)
    tx = (px - cxs[cellx]) / (cxs[cellx + 1L] - cxs[cellx])
    ty = (py - cys[celly]) / (cys[celly + 1L] - cys[celly])
    li = cellx + (celly - 1L) * nc

    pred = function(R){
      (1 - tx) * (1 - ty) * R[li] + tx * (1 - ty) * R[li + 1L] +
        (1 - tx) * ty * R[li + nc] + tx * ty * R[li + nc + 1L]
    }
    maxerr = max(abs(fx$b[1] + fx$b[2] * px + fx$b[3] * py + pred(fx$r) - po[, 1]),
                 abs(fy$b[1] + fy$b[2] * px + fy$b[3] * py + pred(fy$r) - po[, 2]))

    if(is.finite(maxerr) && maxerr < tol){
      F = .warpfield_interp_cpp(fx$r, fy$r, as.integer(cxs), as.integer(cys),
                                fx$b[1], fx$b[2], fx$b[3],
                                fy$b[1], fy$b[2], fy$b[3],
                                as.integer(nx), as.integer(ny))
      # interp returns one nx*ny*1*2 buffer with both planes adjacent, so this
      # is the same cimg imappend(as.cimg(...), as.cimg(...), 'c') would give,
      # without the extra copy.
      fld = imager::as.cimg(F)
      return(list(warpfield = fld, step = s, maxerr = maxerr))
    }

    if(s <= 2L) break
    s = max(2L, s %/% 2L)
  }

  NULL
}

propaneWarp = function(image_in, keyvalues_out=NULL, keyvalues_in=NULL, dim_out = NULL,
                       direction = "auto", boundary = "dirichlet", interpolation = "cubic",
                       doscale = TRUE, dofinenorm = TRUE, plot = FALSE, dotightcrop = TRUE,
                       keepcrop = FALSE, extratight = FALSE, WCSref_out = NULL, WCSref_in = NULL,
                       magzero_out = NULL, magzero_in = NULL, blank = NA, warpfield = NULL,
                       warpfield_return = FALSE, cores = 1, checkWCSequal = FALSE,
                       warpgrid = 'exact', warptol = 1e-5, ...)
{
  if(!requireNamespace("Rwcs", quietly = TRUE)){
    stop("The Rwcs package is needed for this function to work. Please install it from GitHub asgr/Rwcs", call. = FALSE)
  }

  if (!requireNamespace("imager", quietly = TRUE)) {
    stop("The imager package is needed for this function to work. Please install it from CRAN.", call. = FALSE)
  }

  if(!requireNamespace("Rfits", quietly = TRUE)){
    stop("The Rfits package is needed for this function to work. Please install it from GitHub asgr/Rfits", call. = FALSE)
  }

  if(! inherits(image_in, c('Rfits_image', 'Rfits_pointer', 'matrix'))){
    stop('image_in must be either Rfits_image, Rfits_pointer or matrix!')
  }

  if(is.matrix(image_in) & is.null(keyvalues_in)){
    stop('If image_in is a matrix then keyvalues_in must be provided!')
  }

  if(!is.null(keyvalues_in) & is.matrix(image_in)){
    keyvalues_in = keyvalues_in[!is.na(keyvalues_in)]
    image_in = Rfits_create_image(image_in, keyvalues_in)
  }

  keyvalues_in = image_in$keyvalues

  if(any(is.na(keyvalues_in))){
    keyvalues_in = keyvalues_in[!is.na(keyvalues_in)]
    class(keyvalues_in) = 'Rfits_keylist'
  }

  header_in = Rfits_keyvalues_to_raw(keyvalues_in)

  if(any(is.na(keyvalues_out))){
    keyvalues_out = keyvalues_out[!is.na(keyvalues_out)]
    class(keyvalues_out) = 'Rfits_keylist'
  }

  header_out = Rfits_header_to_raw(Rfits_keyvalues_to_header(keyvalues_out))

  if(checkWCSequal){
    if(Rfits_key_match(keyvalues_out, keyvalues_in,
                       check = c('NAXIS1',
                                 'NAXIS2',
                                 'CRPIX1',
                                 'CRPIX2',
                                 'CRVAL1',
                                 'CRVAL2',
                                 'CTYPE1',
                                 'CTYPE2',
                                 'CUNIT1',
                                 'CUNIT1',
                                 'CD1_1',
                                 'CD1_2',
                                 'CD2_1',
                                 'CD2_2'
                        )
                      )
    ){
      message('WCS appears to be the same, directly returning input!')
      image_in$keyvalues$XCUTLO = 1L
      image_in$keyvalues$XCUTHI = dim(image_in)[1]
      image_in$keyvalues$YCUTLO = 1L
      image_in$keyvalues$YCUTHI = dim(image_in)[2]

      image_in$hdr = Rfits_keyvalues_to_hdr(image_in$keyvalues)
      image_in$header = Rfits_keyvalues_to_header(image_in$keyvalues)
      image_in$raw = Rfits_header_to_raw(Rfits_keyvalues_to_header(image_in$keyvalues))

      image_in$keynames = names(image_in$keyvalues)

      image_in$keycomments$XCUTLO = 'Low image x range'
      image_in$keycomments$XCUTHI = 'High image x range'
      image_in$keycomments$YCUTLO = 'Low image y range'
      image_in$keycomments$YCUTHI = 'High image y range'

      return(invisible(image_in))
    }
  }

  if(!is.null(keyvalues_out) & is.null(dim_out)){
    if(isTRUE(keyvalues_out$ZIMAGE)){
      NAXIS1 = keyvalues_out$ZNAXIS1
      NAXIS2 = keyvalues_out$ZNAXIS2
    }else{
      NAXIS1 = keyvalues_out$NAXIS1
      NAXIS2 = keyvalues_out$NAXIS2
    }
    dim_out = c(NAXIS1, NAXIS2)
  }else{
    keyvalues_out$NAXIS1 = dim_out[1]
    keyvalues_out$NAXIS2 = dim_out[2]
  }

  if(is.null(dim_out)){
    stop('Missing NAXIS1 / NAXIS2 in header keyvalues! Specify dim_out.')
  }

  if(dotightcrop){
    suppressMessages({
      BL_out = Rwcs_p2s(0, 0, header=header_out, pixcen='R', WCSref=WCSref_out)
      TL_out = Rwcs_p2s(0, dim_out[2], header=header_out, pixcen='R', WCSref=WCSref_out)
      TR_out = Rwcs_p2s(dim_out[1], dim_out[2], header=header_out, pixcen='R', WCSref=WCSref_out)
      BR_out = Rwcs_p2s(dim_out[1], 0, header=header_out, pixcen='R', WCSref=WCSref_out)
    })

    corners_out = rbind(BL_out, TL_out, TR_out, BR_out)

    suppressMessages({
      tightcrop_out = ceiling(Rwcs_s2p(corners_out, header=header_in, pixcen='R', WCSref=WCSref_in))
    })

    if(is.na(tightcrop_out[1,1])){tightcrop_out[1,] = c(0,0)}
    if(is.na(tightcrop_out[2,1])){tightcrop_out[2,] = c(0,dim(image_in)[2])}
    if(is.na(tightcrop_out[3,1])){tightcrop_out[3,] = c(dim(image_in)[1],dim(image_in)[2])}
    if(is.na(tightcrop_out[4,1])){tightcrop_out[4,] = c(dim(image_in)[1],0)}

    min_x_out = max(1L, min(tightcrop_out[,1]))
    max_x_out = min(dim(image_in)[1], max(tightcrop_out[,1]))
    min_y_out = max(1L, min(tightcrop_out[,2]))
    max_y_out = min(dim(image_in)[2], max(tightcrop_out[,2]))

    # An inverted range means the output frame barely (or not at all) overlaps
    # the input. Handing c(hi, lo) to the crop below is a *box* crop, so it
    # silently expands the working grid rather than failing -- measured up to
    # 4.5x the pixel count, with a normal-sized but all-NaN result and no
    # message. Warn, but otherwise leave the arithmetic untouched.
    if (min_x_out > max_x_out || min_y_out > max_y_out) {
      warning(sprintf(
        paste0('tight crop is inverted (x [%d, %d], y [%d, %d] over a %dx%d input): ',
               'the output frame has little or no overlap with the input. ',
               'Warped data is expected to be empty.'),
        min_x_out, max_x_out, min_y_out, max_y_out,
        dim(image_in)[1], dim(image_in)[2]),
        call. = FALSE)
    }

    if(min_x_out != 1 | max_x_out != dim(image_in)[1] | min_y_out != 1 | max_y_out != dim(image_in)[2]){
      if(inherits(image_in, 'Rfits_pointer')){
        image_in = image_in[c(min_x_out, max_x_out), c(min_y_out, max_y_out), header=TRUE]
      }else{
        image_in = image_in[c(min_x_out, max_x_out), c(min_y_out, max_y_out)]
      }

      keyvalues_in = image_in$keyvalues
      header_in = Rfits_header_to_raw(Rfits_keyvalues_to_header(keyvalues_in))
    }else{
      if(inherits(image_in, 'Rfits_pointer')){
        image_in = image_in[,]
      }
    }

    suppressMessages({
      BL_in = Rwcs_p2s(0, 0,header=header_in, pixcen='R', WCSref=WCSref_in)
      TL_in = Rwcs_p2s(0, dim(image_in)[2], header=header_in, pixcen='R', WCSref=WCSref_in)
      TR_in = Rwcs_p2s(dim(image_in)[1], dim(image_in)[2], header=header_in, pixcen='R', WCSref=WCSref_in)
      BR_in = Rwcs_p2s(dim(image_in)[1], 0, header=header_in, pixcen='R', WCSref=WCSref_in)
    })

    corners_in = rbind(BL_in, TL_in, TR_in, BR_in)

    suppressMessages({
      tightcrop_in = ceiling(Rwcs_s2p(corners_in, header=header_out, pixcen='R', WCSref=WCSref_out))
    })

    if(is.na(tightcrop_in[1,1])){tightcrop_in[1,] = c(0,0)}
    if(is.na(tightcrop_in[2,1])){tightcrop_in[2,] = c(0,dim_out[2])}
    if(is.na(tightcrop_in[3,1])){tightcrop_in[3,] = c(dim_out[1],dim_out[2])}
    if(is.na(tightcrop_in[4,1])){tightcrop_in[4,] = c(dim_out[1],0)}

    min_x_in = max(1L, min(tightcrop_in[,1]))
    max_x_in = max(min_x_in + dim(image_in)[1] - 1L, range(tightcrop_in[,1])[2])
    min_y_in = max(1L, min(tightcrop_in[,2]))
    max_y_in = max(min_y_in + dim(image_in)[2] - 1L, range(tightcrop_in[,2])[2])

    # new code should be more efficient!

    if(isTRUE(keyvalues_out$ZIMAGE)){
      keyvalues_out$ZNAXIS1 = max_x_in - min_x_in + 1L
      keyvalues_out$ZNAXIS2 = max_y_in - min_y_in + 1L
    }else{
      keyvalues_out$NAXIS1 = max_x_in - min_x_in + 1L
      keyvalues_out$NAXIS2 = max_y_in - min_y_in + 1L
    }

    keyvalues_out$CRPIX1 = keyvalues_out$CRPIX1 - min_x_in + 1L
    keyvalues_out$CRPIX2 = keyvalues_out$CRPIX2 - min_y_in + 1L

    header_out = Rfits_header_to_raw(Rfits_keyvalues_to_header(keyvalues_out))

    image_out = list(
      imDat = matrix(c(blank,image_in$imDat[0]), max_x_in - min_x_in + 1L, max_y_in - min_y_in + 1L),
      keyvalues = keyvalues_out,
      hdr = Rfits_keyvalues_to_hdr(keyvalues_out),
      header = Rfits_keyvalues_to_header(keyvalues_out),
      raw = Rfits_header_to_raw(Rfits_keyvalues_to_header(keyvalues_out)),
      keynames = names(keyvalues_out),
      keycomments = as.list(rep('', length(keyvalues_out)))
    )
    names(image_out$keycomments) = image_out$keynames
    class(image_out) = c('Rfits_image', class(image_out))
  }else{
    if(inherits(image_in, 'Rfits_pointer')){
      image_in = image_in[,]
    }

    image_out = list(
      imDat = matrix(c(blank,image_in$imDat[0]), max(dim(image_in)[1], dim_out[1]), max(dim(image_in)[2], dim_out[2])),
      keyvalues = keyvalues_out,
      hdr = Rfits_keyvalues_to_hdr(keyvalues_out),
      header = Rfits_keyvalues_to_header(keyvalues_out),
      raw = Rfits_header_to_raw(Rfits_keyvalues_to_header(keyvalues_out)),
      keynames = names(keyvalues_out),
      keycomments = as.list(rep('', length(keyvalues_out)))
    )
    names(image_out$keycomments) = image_out$keynames
    class(image_out) = c('Rfits_image', class(image_out))

    min_x_in = 1L
    max_x_in = dim_out[1]
    min_y_in = 1L
    max_y_in = dim_out[2]
  }

  dim_min_x_in = min(dim(image_in)[1], dim(image_out$imDat)[1])
  dim_min_y_in = min(dim(image_in)[2], dim(image_out$imDat)[2])

  if(anyInfinite(image_in$imDat)){
    image_in$imDat[is.infinite(image_in$imDat)] = NA
  }

  image_out$imDat[1:dim_min_x_in, 1:dim_min_y_in] = image_in$imDat[1:dim_min_x_in, 1:dim_min_y_in]
  rm(image_in)

  if(!is.null(magzero_in) & !is.null(magzero_out)){
    image_out$imDat = image_out$imDat*10^(-0.4*(magzero_in - magzero_out))
    keyvalues_out$MAGZERO = magzero_out
  }

  suppressMessages({
    pixscale_in = pixscale(keyvalues_in)
    pixscale_out = pixscale(keyvalues_out)
  })

  if (direction == "auto") {
    if (pixscale_in < pixscale_out) {
      direction = "forward"
    }
    if (pixscale_in >= pixscale_out) {
      direction = "backward"
    }
  }

  if(is.null(warpfield)){
    dim_field = dim(image_out$imDat)[1:2]

    warpfun = if (direction == "forward") {
      function(x, y, cores, ...) .warpfunc_in2out(
        x = x, y = y,
        header_in = header_in, WCSref_in = WCSref_in,
        header_out = header_out, WCSref_out = WCSref_out,
        cores = cores)
    } else {
      function(x, y, cores, ...) .warpfunc_out2in(
        x = x, y = y,
        header_in = header_in, WCSref_in = WCSref_in,
        header_out = header_out, WCSref_out = WCSref_out,
        cores = cores)
    }

    built = NULL
    if(!identical(warpgrid, 'exact')){
      step0 = if(is.numeric(warpgrid)) as.integer(warpgrid[1]) else 64L
      built = .warpfield_coarse(warpfun, dim_field, tol = warptol,
                                step0 = step0, cores = cores)
      if(is.null(built)){
        message('coarse warpgrid did not converge to tolerance; using exact field.')
      }else{
        message(sprintf('coarse warpgrid: step %d (%d pts), max field error %.2e px',
                        built$step,
                        length(seq(1, dim_field[1], by = built$step)) *
                          length(seq(1, dim_field[2], by = built$step)),
                        built$maxerr))
        warpfield = built$warpfield
      }
    }

    if(is.null(built)){
      pix_grid = expand.grid(1:dim_field[1], 1:dim_field[2])
      warp_out = warpfun(pix_grid[, 1], pix_grid[, 2], cores = cores)

      warpmat1 = matrix(warp_out[, 1], dim_field[1], dim_field[2])

      if(anyInfinite(warpmat1)){
        message('Infinity found in warpfield- patching!')
        warpmat1[is.infinite(warpmat1)] = NA
        warpmat1 = propanePatchPix(warpmat1)
      }

      warpmat2 = matrix(warp_out[, 2], dim_field[1], dim_field[2])

      if(anyInfinite(warpmat2)){
        message('Infinity found in warpfield- patching!')
        warpmat2[is.infinite(warpmat2)] = NA
        warpmat2 = propanePatchPix(warpmat2)
      }

      warpfield = imager::imappend(list(
        imager::as.cimg(warpmat1),
        imager::as.cimg(warpmat2)
      ), 'c')

      rm(pix_grid)
      rm(warp_out)
      rm(warpmat1)
      rm(warpmat2)
    }
  }

  image_out$imDat = imager::warp(
    im = imager::as.cimg(image_out$imDat),
    warpfield = warpfield,
    mode = switch(direction, backward = 0L, forward =
                    2L),
    interpolation = switch(
      interpolation,
      nearest = 0L,
      linear = 1L,
      cubic = 2L
    ),
    boundary_conditions = switch(
      boundary,
      dirichlet = 0L,
      neumann = 1L,
      periodic = 2L
    )
  )

  if (dofinenorm) {
    norm = matrix(1, dim(image_out$imDat)[1], dim(image_out$imDat)[2])
    norm = imager::warp(
      im = imager::as.cimg(norm),
      warpfield = warpfield,
      mode = switch(direction, backward = 0L, forward = 2L),
      interpolation = switch(
        interpolation,
        nearest = 0L,
        linear = 1L,
        cubic = 2L
      ),
      boundary_conditions = switch(
        boundary,
        dirichlet = 0L,
        neumann = 1L,
        periodic = 2L
      )
    )

    image_out$imDat = image_out$imDat / norm
    rm(norm)
  }

  if (doscale) {
    image_out$imDat = image_out$imDat * (pixscale_out / pixscale_in) ^ 2
  }

  image_out$imDat = as.matrix(image_out$imDat)

  if(dotightcrop==FALSE | keepcrop==FALSE){
    image_out = image_out[c(1L - (min_x_in - 1L), dim_out[1] - (min_x_in - 1L)),c(1L - (min_y_in - 1L), dim_out[2] - (min_y_in - 1L)), box=1] #box=1 just in case we have a single pixel left
    image_out$keyvalues$XCUTLO = 1L
    image_out$keyvalues$XCUTHI = dim_out[1]
    image_out$keyvalues$YCUTLO = 1L
    image_out$keyvalues$YCUTHI = dim_out[2]

    image_out$hdr = Rfits_keyvalues_to_hdr(image_out$keyvalues)
    image_out$header = Rfits_keyvalues_to_header(image_out$keyvalues)
    image_out$raw = Rfits_header_to_raw(Rfits_keyvalues_to_header(image_out$keyvalues))

    image_out$keynames = names(image_out$keyvalues)

    image_out$keycomments$XCUTLO = 'Low image x range'
    image_out$keycomments$XCUTHI = 'High image x range'
    image_out$keycomments$YCUTLO = 'Low image y range'
    image_out$keycomments$YCUTHI = 'High image y range'
  }else{

    if(max_x_in > dim_out[1]){
      trim_x = max_x_in - dim_out[1]
      if(dim(image_out)[1] > trim_x){
        image_out = image_out[1:(dim(image_out)[1] - trim_x), , box=1] #box=1 just in case we have a single pixel left
      }else{
        image_out = image_out[1, , box=1] #box=1 just in case we have a single pixel left
      }
      max_x_in = dim_out[1]
    }

    if(max_y_in > dim_out[2]){
      trim_y = max_y_in - dim_out[2]
      if(dim(image_out)[2] > trim_y){
        image_out = image_out[, 1:(dim(image_out)[2] - trim_y), box=1] #box=1 just in case we have a single pixel left
      }else{
        image_out = image_out[, 1, box=1] #box=1 just in case we have a single pixel left
      }
      max_y_in = dim_out[2]
    }

    if(extratight){
      final_pix = which(!is.na(image_out$imDat), arr.ind = TRUE)
      crop_x_lo = min(final_pix[,1])
      crop_x_hi = max(final_pix[,1])
      crop_y_lo = min(final_pix[,2])
      crop_y_hi = max(final_pix[,2])

      #final extra tight crop
      image_out = image_out[c(crop_x_lo, crop_x_hi), c(crop_y_lo, crop_y_hi)]

      #update where we are in the parent frame
      min_x_in = min_x_in + crop_x_lo - 1L
      max_x_in = min_x_in + dim(image_out)[1] -1L
      min_y_in = min_y_in + crop_y_lo - 1L
      max_y_in = min_y_in + dim(image_out)[2] -1L
    }

    image_out$keyvalues$XCUTLO = min_x_in
    image_out$keyvalues$XCUTHI = max_x_in
    image_out$keyvalues$YCUTLO = min_y_in
    image_out$keyvalues$YCUTHI = max_y_in

    image_out$keycomments$XCUTLO = 'Low image x range'
    image_out$keycomments$XCUTHI = 'High image x range'
    image_out$keycomments$YCUTLO = 'Low image y range'
    image_out$keycomments$YCUTHI = 'High image y range'

    image_out$keynames = names(image_out$keyvalues)

    image_out$crop = c(xlo=min_x_in, xhi=max_x_in, ylo=min_y_in, yhi=max_y_in) #we want to keep the subset location for potential later writing
  }

  image_out$history = c(image_out$history, "Warped with propaneWarp")

  image_out = Rfits_check_image(image_out)

  if(plot){
    plot(image_out, ...)
  }

  if(warpfield_return){
    image_out$warpfield = warpfield
  }

  return(invisible(image_out))
}

propaneRebin = function(image, scale = 1,interpolation = 6){
  if (!requireNamespace("imager", quietly = TRUE)) {
    stop("The imager package is needed for this function to work. Please install it from CRAN.",
         call. = FALSE)
  }

  imdim = dim(image)
  if(scale > 1){
    #scale = floor(scale)
    size_x = floor(imdim[1]*scale - (scale - 1L))
    size_y = floor(imdim[2]*scale - (scale - 1L))
  }else{
    #scale = 1/ceiling(1/scale)
    size_x = floor(imdim[1]*scale)
    size_y = floor(imdim[2]*scale)
  }

  if(inherits(image,'Rfits_image')){
    image_resize = as.matrix(imager::resize(im=imager::as.cimg(image$imDat), size_x=size_x, size_y=size_y, interpolation_type=interpolation))
    norm = matrix(1, dim(image$imDat)[1], dim(image$imDat)[2])
    norm_resize = as.matrix(imager::resize(im=imager::as.cimg(norm), size_x=size_x, size_y=size_y, interpolation_type=interpolation))
    image_resize = (image_resize / norm_resize) / scale^2

    keyvalues_out = image$keyvalues

    if(isTRUE(image$keyvalues$ZIMAGE)){
      keyvalues_out$ZNAXIS1 = dim(image_resize)[1]
      keyvalues_out$ZNAXIS2 = dim(image_resize)[2]
    }else{
      keyvalues_out$NAXIS1 = dim(image_resize)[1]
      keyvalues_out$NAXIS2 = dim(image_resize)[2]
    }

    if(scale > 1){
      keyvalues_out$CRPIX1 = (keyvalues_out$CRPIX1 - 0.5) * scale
      keyvalues_out$CRPIX2 = (keyvalues_out$CRPIX2 - 0.5) * scale
    }else{
      keyvalues_out$CRPIX1 = (keyvalues_out$CRPIX1 + 0.5) * scale
      keyvalues_out$CRPIX2 = (keyvalues_out$CRPIX2 + 0.5) * scale
    }
    keyvalues_out$CD1_1 = keyvalues_out$CD1_1 / scale
    keyvalues_out$CD1_2 = keyvalues_out$CD1_2 / scale
    keyvalues_out$CD2_1 = keyvalues_out$CD2_1 / scale
    keyvalues_out$CD2_2 = keyvalues_out$CD2_2 / scale

    image_out = list(
      imDat = image_resize,
      keyvalues = keyvalues_out,
      hdr = Rfits_keyvalues_to_hdr(keyvalues_out),
      header = Rfits_keyvalues_to_header(keyvalues_out),
      raw = Rfits_header_to_raw(Rfits_keyvalues_to_header(keyvalues_out)),
      keynames = names(keyvalues_out),
      keycomments = as.list(rep('', length(keyvalues_out)))
    )
    names(image_out$keycomments) = image_out$keynames
    class(image_out) = c('Rfits_image', class(image_out))
    image_out = Rfits_check_image(image_out)

    return(image_out)
  }else{
    image_resize = as.matrix(imager::resize(im=imager::as.cimg(image), size_x=size_x, size_y=size_y, interpolation_type=interpolation))
    norm = matrix(1, dim(image)[1], dim(image)[2])
    norm_resize = as.matrix(imager::resize(im=imager::as.cimg(norm), size_x=size_x, size_y=size_y, interpolation_type=interpolation))
    image_resize = (image_resize / norm_resize) / scale^2

    return(image_resize)
  }
}

# Can a warp field built for one band be reused for another? The field depends
# only on the input pixel grid, the input WCS, the output WCS and the direction
# -- never on the band's pixel values -- so bands of one detection can share a
# single field. Only the WCS-bearing keys matter here; the full keyvalues lists
# differ between bands in ways that are irrelevant to the geometry (EXTNAME,
# MAGZERO, filter, ...) and must not defeat the comparison.
WCS_GEOM_PATTERN = paste0(
  '^(CTYPE|CRVAL|CRPIX|CD[12]_|PC[12]_|CDELT|CUNIT|CROTA|LONPOLE|LATPOLE|',
  'RADE|EQUINOX|RADESYS|WCSAXES|NAXIS|PV[12]_|A_[12]|B_[12]|D_[12]|',
  'SIP|POLORDER|ZP[12]|ZIMAGE|ZNAXIS|ZCRPIX|ZCD[12]_)')

.warpfield_geom = function(image){
  kv = image$keyvalues
  sel = grep(WCS_GEOM_PATTERN, names(kv), value = TRUE)
  list(dim = dim(image)[1:2],
       geom = paste(names(kv)[names(kv) %in% sel], unname(kv)[names(kv) %in% sel],
                    sep = '=', collapse = '\n'))
}

.warpfield_reusable = function(field, geom_prev, image_next){
  if(is.null(field) || !inherits(field, 'cimg')) return(FALSE)
  if(is.null(geom_prev)) return(FALSE)
  geom_next = .warpfield_geom(image_next)
  identical(geom_next$dim, geom_prev$dim) &&
    identical(geom_next$geom, geom_prev$geom)
}

propaneWarpProPane = function(propane_in, keyvalues_out=NULL, dim_out = NULL, magzero_out = NULL, ..., warpfield_share = TRUE){

  if(!is.null(magzero_out)){
    zero_point_scale = 10^(-0.4*(propane_in$image$keyvalues$MAGZERO - magzero_out))
  }else{
    magzero_out = propane_in$image$keyvalues$MAGZERO
    zero_point_scale = 1
  }

  keyvalues_out$R_VER = propane_in$image$keyvalues$R_VER
  keyvalues_out$PANE_VER = propane_in$image$keyvalues$PANE_VER
  keyvalues_out$RWCS_VER = propane_in$image$keyvalues$RWCS_VER

  # The warp field depends only on the input grid, the input/output WCS and the
  # direction -- none of which vary between bands of the same detection. Build
  # it once and hand it to the rest. Passing warpfield= explicitly in ... opts
  # out of the sharing logic entirely.
  dots = list(...)
  shared_field = NULL
  shared_geom = NULL
  # An explicit field, or a request to see the fields, opts out of sharing.
  manual = !is.null(dots$warpfield) || isTRUE(dots$warpfield_return)

  warp_band = function(band, doscale){
    args = c(list(image_in = band, keyvalues_out = keyvalues_out,
                  dim_out = dim_out, doscale = doscale), dots)

    share = isTRUE(warpfield_share) && !manual
    if(share){
      if(.warpfield_reusable(shared_field, shared_geom, band)){
        args$warpfield = shared_field
      }else{
        # First usable band, or one whose geometry differs: build on this call.
        args$warpfield_return = TRUE
      }
    }

    out = do.call(propaneWarp, args)

    if(share && !is.null(out$warpfield)){
      if(is.null(shared_field)){
        shared_field <<- out$warpfield
        shared_geom <<- .warpfield_geom(band)
      }
      # Keep the returned structure identical to the unshared path.
      out$warpfield = NULL
    }
    out
  }

  if(!is.null(propane_in$image)){
    message('warping image')

    image_warp = warp_band(propane_in$image*zero_point_scale, doscale = TRUE)

    image_warp$keyvalues$EXTNAME = 'image'
    image_warp$keyvalues$MAGZERO = magzero_out
    image_warp$history = c(propane_in$image$history, image_warp$history)
    image_warp = Rfits_check_image(image_warp)
  }else{
    image_warp = NULL
  }

  if(!is.null(propane_in$weight)){
    message('warping weight')

    weight_warp = warp_band(propane_in$weight, doscale = FALSE)

    weight_warp$keyvalues$EXTNAME = 'weight'
    weight_warp$history = c(propane_in$weight$history, weight_warp$history)
    weight_warp = Rfits_check_image(weight_warp)
  }else{
    weight_warp = NULL
  }

  if(!is.null(propane_in$inVar)){
    message('warping inVar')

    inVar_warp = warp_band(propane_in$inVar/(zero_point_scale^2), doscale = FALSE)*
      (pixscale(propane_in$inVar$keyvalues)^4 / pixscale(keyvalues_out)^4)

    inVar_warp$keyvalues$EXTNAME = 'inVar'
    inVar_warp$keyvalues$MAGZERO = magzero_out
    inVar_warp$history = c(propane_in$inVar$history, inVar_warp$history)
    inVar_warp = Rfits_check_image(inVar_warp)
  }else{
    inVar_warp = NULL
  }

  if(!is.null(propane_in$exp)){
    message('warping exp')

    exp_warp = warp_band(propane_in$exp, doscale = FALSE)

    exp_warp$keyvalues$EXTNAME = 'exp'
    exp_warp$history = c(propane_in$exp$history, exp_warp$history)
    exp_warp = Rfits_check_image(exp_warp)
  }else{
    exp_warp = NULL
  }

  if(!is.null(propane_in$cold)){
    message('warping cold')

    cold_warp = warp_band(propane_in$cold*zero_point_scale, doscale = TRUE)

    cold_warp$keyvalues$EXTNAME = 'cold'
    cold_warp$keyvalues$MAGZERO = magzero_out
    cold_warp$history = c(propane_in$cold$history, cold_warp$history)
    cold_warp = Rfits_check_image(cold_warp)
  }else{
    cold_warp = NULL
  }

  if(!is.null(propane_in$hot)){
    message('warping hot')

    hot_warp = warp_band(propane_in$hot*zero_point_scale, doscale = TRUE)

    hot_warp$keyvalues$EXTNAME = 'hot'
    hot_warp$keyvalues$MAGZERO = magzero_out
    hot_warp$history = c(propane_in$hot$history, hot_warp$history)
    hot_warp = Rfits_check_image(hot_warp)
  }else{
    hot_warp = NULL
  }

  if(!is.null(propane_in$clip)){
    message('warping clip')

    clip_warp = warp_band(propane_in$clip, doscale = FALSE)

    clip_warp$keyvalues$EXTNAME = 'clip'
    clip_warp$history = c(propane_in$clip$history, clip_warp$history)
    clip_warp = Rfits_check_image(clip_warp)
  }else{
    clip_warp = NULL
  }

  output = list(
    image = image_warp,
    weight = weight_warp,
    inVar = inVar_warp,
    exp = exp_warp,
    cold = cold_warp,
    hot = hot_warp,
    clip = clip_warp
  )

  class(output) = "ProPane"
  return(invisible(output))
}

propaneSegimWarp = function(segim=NULL, ...){
  invisible(propaneWarp(image_in=segim, direction='backward', interpolation='nearest', doscale=FALSE, ...)$imDat)
}

propaneSegimShare=function(segim=NULL, keyvalues_in=NULL, keyvalues_out=NULL, pixcut=1){
  segID_in = sort(unique(as.integer(segim[segim > 0])))
  segim_warp = propaneSegimWarp(segim = segim,
                                 keyvalues_in = keyvalues_in,
                                 keyvalues_out = keyvalues_out)
  segim_warp_tab = tabulate(segim_warp)
  segID_warp = which(segim_warp_tab >= pixcut)
  segim_warp[!segim_warp %in% segID_warp] = 0
  segim_unwarp = propaneSegimWarp(segim = segim_warp,
                                   keyvalues_in = keyvalues_out,
                                   keyvalues_out = keyvalues_in)
  segim_out = NULL
  segimDT = data.table(segim = as.integer(segim),
                       segim_out = as.integer(segim_unwarp))
  segimDT = segimDT[segim > 0, ]
  segim_groups = segimDT[, list(segim_out = list(tabulate(segim_out))), keyby =
                           segim]
  sharemat = matrix(0, max(segim_warp), dim(segim_groups)[1])
  for (i in 1:(dim(sharemat)[2])) {
    sharemat[1:length(unlist(segim_groups$segim_out[i])), i] = unlist(segim_groups$segim_out[i])
  }
  sharemat = sharemat / rowSums(sharemat)
  sharemat = sharemat[is.finite(sharemat[, 1]), , drop = FALSE]
  colnames(sharemat) = segID_in
  rownames(sharemat) = segID_warp
  # if (!is.null(weights)) {
  #   t(t(sharemat) * weights)
  #   sharemat = sharemat / rowSums(sharemat)
  # }
  shareseg = diag(sharemat[segID_warp %in% segID_in, segID_in %in% segID_warp])
  invisible(
    list(
      segID_in = segID_in,
      segID_warp = segID_warp,
      segim_warp = segim_warp,
      sharemat = sharemat,
      shareseg = shareseg
    )
  )
}
