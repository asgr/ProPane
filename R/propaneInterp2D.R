propaneInterp2D = function(x, y, z = 1, image, xlim=NULL, ylim=NULL, pixcen='R', type = 'add', zero = FALSE) {
  if(missing(x)){
    stop("Must provide x")
  }
  if(length(dim(x) == 2)){
    if(dim(x)[2] == 2){
      y = x[,2]
      x = x[,1]
    }else if(dim(x)[2] == 3){
      z = x[,3]
      y = x[,2]
      x = x[,1]
    }else{
      stop("x must only have 2 or 3 columns")
    }
  }
  if(length(x) != length(y)){
    stop('x and y must be same length')
  }
  if(length(x) != length(z) & length(z) != 1){
    stop('z and x/y must be same length, or length 1')
  }
  if(missing(image)){
    stop("Need image matrix")
  }

  if(!is.null(xlim)){
    xseq = seq(xlim[1], xlim[2], len=dim(image)[1])
    x = (dim(image)[1] - 1)*(x - xlim[1])/diff(range(xlim))
    if(pixcen == 'FITS'){
      xseq = xseq - 0.5
      x = x + 1
    }else{
      x = x + 0.5
    }
  }else{
    if(!is.null(ylim)){
      xseq = 1:dim(image)[1] - 0.5
    }
    xseq = NULL
  }

  if(!is.null(ylim)){
    yseq = seq(ylim[1], ylim[2], len=dim(image)[2])
    y = (dim(image)[2] - 1)*(y - ylim[1])/diff(range(ylim))
    if(pixcen == 'FITS'){
      y = y + 1
      yseq = yseq - 0.5
    }else{
      y = y + 0.5
    }
  }else{
    if(!is.null(xlim)){
      yseq = 1:dim(image)[2] - 0.5
    }
    yseq = NULL
  }

  image = .propaneInterp2D(
                   x = x,
                   y = y,
                   z = z,
                   image = image,
                   FITS = switch(pixcen, FITS = TRUE, R = FALSE),
                   type = switch(type, add = 1L, sub = 2L),
                   zero = zero
                   )
  if(!is.null(xseq) & !is.null(yseq)){
    image = list(x=xseq, y=yseq, z=image)
  }

  return(image)
}

propaneBin2D = function(x, y, z = 1, image, xlim=NULL, ylim=NULL, pixcen='R', type = 'add', zero = FALSE) {
  if(missing(x)){
    stop("Must provide x")
  }
  if(length(dim(x) == 2)){
    if(dim(x)[2] == 2){
      y = x[,2]
      x = x[,1]
    }else if(dim(x)[2] == 3){
      z = x[,3]
      y = x[,2]
      x = x[,1]
    }else{
      stop("x must only have 2 or 3 columns")
    }
  }
  if(length(x) != length(y)){
    stop('x and y must be same length')
  }
  if(length(x) != length(z) & length(z) != 1){
    stop('z and x/y must be same length, or length 1')
  }
  if(missing(image)){
    stop("Need image matrix")
  }

  if(!is.null(xlim)){
    xseq = seq(xlim[1], xlim[2], len=dim(image)[1])
    x = (dim(image)[1] - 1)*(x - xlim[1])/diff(range(xlim))
    if(pixcen == 'FITS'){
      xseq = xseq - 0.5*diff(range(xlim))/dim(image)[1]
      x = x + 1
    }else{
      x = x + 0.5
    }
  }else{
    if(!is.null(ylim)){
      xseq = 1:dim(image)[1] - 0.5
    }else{
      xseq = NULL
    }
  }

  if(!is.null(ylim)){
    yseq = seq(ylim[1], ylim[2], len=dim(image)[2])
    y = (dim(image)[2] - 1)*(y - ylim[1])/diff(range(ylim))
    if(pixcen == 'FITS'){
      y = y + 1
      yseq = yseq - 0.5*diff(range(ylim))/dim(image)[2]
    }else{
      y = y + 0.5
    }
  }else{
    if(!is.null(xlim)){
      yseq = 1:dim(image)[2] - 0.5
    }else{
      yseq = NULL
    }
  }

  if(length(z) == 1 & z[1] == 1){
    .propaneBin2Dint(
      x = x,
      y = y,
      image = image,
      FITS = switch(pixcen, FITS = TRUE, R = FALSE),
      type = switch(type, add = 1L, sub = 2L),
      zero = zero
    )
  }else{
    .propaneBin2D(
      x = x,
      y = y,
      z = z,
      image = image,
      FITS = switch(pixcen, FITS = TRUE, R = FALSE),
      type = switch(type, add = 1L, sub = 2L),
      zero = zero
    )
  }

  if(!is.null(xseq) & !is.null(yseq)){
    image = list(x=xseq, y=yseq, z=image)
  }

  return(image)
}
