propaneLocalFunc = function(image, imager_func=NULL, dither=1, offset=1, iter=1, kern='square',
                            cores=1, multitype='fork', verbose=TRUE, ...){

  if(missing(image)){
    stop('Need image!')
  }

  if(!is.matrix(image)){
    stop('image must be a matrix!')
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

    image_list = foreach(i = 1:dim(cutgrid)[1])%dopar%{
      temp = magcutout(
                image = image,
                loc = c(dim(image)[1]/2 + cutgrid[i,1], dim(image)[2]/2 + cutgrid[i,2]),
                box = dim(image)
                )$image
      return(imager::as.cimg(temp))
    }

    image = propaneStackFlatFunc(image_list, imager_func=imager_func, ...)$image
  }

  return(image)
}
