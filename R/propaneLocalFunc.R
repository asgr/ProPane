propaneLocalMed = function(image, dither=1, iter=1, threshold=Inf, maxdiff=Inf, verbose=TRUE){

  if (!requireNamespace("imager", quietly = TRUE)) {
    stop("The imager package is needed for this function to work. Please install it from CRAN.", call. = FALSE)
  }

  if(missing(image)){
    stop('Need image!')
  }

  if(imager::is.cimg(image)){
    #do nothing
  }else{
    if(is.matrix(image)){
      image = imager::as.cimg(image)
    }else{
      stop('image must be a Cimg or a matrix!')
    }
  }

  if(length(dither) == 1L){
    dither = rep(dither, iter)
  }

  if(length(threshold) == 1L){
    threshold = rep(threshold, iter)
  }

  if(length(maxdiff) == 1L){
    maxdiff = rep(maxdiff, iter)
  }

  for(i in 1:iter){
    if(iter > 1L & verbose){
      message('Iteration ', i,' of ', iter)
    }

    box = dither[i]*2 + 1

    image_new = imager::medianblur(image, n=box, threshold=threshold[i]) #the Inf means NAs will work correctly

    if(is.finite(maxdiff[i])){
      sel = which(image_new - image > maxdiff[i])
      image_new[sel] = image[sel]
    }

    image = image_new
  }

  return(image)
}

propaneLocalFunc = function(image, imager_func=NULL, dither=1, offset=1, iter=1, kern = 'square',
                            cores=1, multitype='fork', ondisk=FALSE, dump_dir=tempdir(),
                            verbose=TRUE, ...){

  if(missing(image)){
    stop('Need image!')
  }

  if(!is.matrix(image)){
    stop('image must be a matrix!')
  }

  if(ondisk){
    dump_dir = path.expand(dump_dir)
    if(dir.exists(dump_dir) == FALSE){
      dir.create(dump_dir, recursive = TRUE)
    }
    message('Frames being dumped to ', dump_dir)
  }

  if(length(dither) == 1L){
    dither = rep(dither, iter)
  }

  if(length(offset) == 1L){
    offset = rep(offset, iter)
  }

  if(multitype=='fork'){
    registerDoParallel(cores=cores)
  }else if(multitype=='cluster'){
    registerDoParallel(cl=cores)
  }

  kern = tolower(kern)

  for(i in 1:iter){
    if(iter > 1L & verbose){
      message('Iteration ', i,' of ', iter)
    }

    boxlim = dither[i]*offset[i]
    cutgrid = expand.grid(seq(-boxlim, boxlim, by=offset[i]), seq(-boxlim, boxlim, by=offset[i]))

    if(kern == 'circle'){
      cutgrid = cutgrid[cutgrid[,1]^2 + cutgrid[,2]^2 <= dither,]
    }else if(kern == 'square'){
      #Do nothing
    }else{
      stop('kern must be one of square of circle!')
    }

    if(verbose){
      message(' - Creating image_list')
    }

    image_list = foreach(i = 1:dim(cutgrid)[1])%dopar%{
      temp = magcutout(
                image = image,
                loc = c(dim(image)[1]/2 + cutgrid[i,1], dim(image)[2]/2 + cutgrid[i,2]),
                box = dim(image)
                )$image

      if(ondisk){
        tempfile = paste0(dump_dir,'/image_local_den_',i,'.fits')
        Rfits_write_image(temp, tempfile)
        return(Rfits_point(tempfile))
      }else{
        return(imager::as.cimg(temp))
      }
    }

    if(verbose){
      message(' - Calculating local median')
    }

    image = propaneStackFlatFunc(image_list,
                                 imager_func = imager_func,
                                 ondisk = ondisk,
                                 cores = cores,
                                 multitype = multitype,
                                 ...)$image
  }

  return(image)
}
