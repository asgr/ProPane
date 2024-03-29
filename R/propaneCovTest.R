propaneCovTest = function(image_in, keyvalues_out=NULL, sparse=7L, flux_map=FALSE, degree=1, ...){
  dim_in = dim(image_in)

  if(sparse == 'auto'){
    sparse = ceiling(pixscale(image_in)/pixscale(keyvalues_out)) #something here!!!!!!
    sparse = max(sparse, 2L)
  }

  x_seq = seq(ceiling(sparse/2), dim_in[1], by=sparse)
  y_seq = seq(ceiling(sparse/2), dim_in[2], by=sparse)

  image_in$imDat[] = 0L

  image_in$imDat[x_seq, y_seq] = 1L

  output = propaneWarp(image_in = image_in,
              keyvalues_out = keyvalues_out,
              doscale = FALSE,
              ...
              )

  hi_frac = 1/sparse^2
  hi_quan = quantile(output$imDat, 1 - hi_frac, na.rm=TRUE)
  hi_mean = mean(output$imDat[output$imDat >= hi_quan], na.rm=TRUE)

  lo_quan = quantile(output$imDat, max(0, 1 - 9*hi_frac), na.rm=TRUE)
  lo_mean = mean(output$imDat[output$imDat >= lo_quan & output$imDat <= hi_quan], na.rm=TRUE)

  if(flux_map){
    im_dim = dim(output$imDat)
    im_pix = expand.grid(x=1:im_dim[1], y=1:im_dim[2])
    sel = which(!is.na(output$imDat))
    im_pix_sub = im_pix[sel,]
    lm_out = lm(output$imDat[sel]*sparse^2 ~ poly(x, degree=degree) + poly(y, degree=degree), data=im_pix_sub)
    flux_map = matrix(predict(lm_out, im_pix), im_dim[1], im_dim[2])
    return(list(cov_im = output, lo_mean=lo_mean, hi_mean=hi_mean, flux_map=flux_map))
  }else{
    return(list(cov_im = output, lo_mean=lo_mean, hi_mean=hi_mean))
  }
}
