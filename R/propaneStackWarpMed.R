propaneStackWarpMed = function(
    filelist = NULL,
    dirlist = NULL,
    extlist = 1,
    pattern = NULL,
    recursive = TRUE,
    zap = NULL,
    keyvalues_out = NULL,
    cores = 4,
    chunk = 1e3,
    doweight = TRUE,
    useCUTLO = TRUE){

  timestart = proc.time()[3]

  j = NULL

  assertList(keyvalues_out, null.ok = TRUE)
  assertIntegerish(cores, len=1)
  assertIntegerish(chunk, len=1)

  if(!requireNamespace("Rfits", quietly = TRUE)){
    stop('The Rfits package is needed for stacking to work. Please install from GitHub asgr/Rfits.', call. = FALSE)
  }

  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for stacking to work. Please install from GitHub asgr/Rfits.', call. = FALSE)
  }

  registerDoParallel(cores=cores)

  image_list = Rfits_make_list(filelist = filelist,
                               dirlist = dirlist,
                               extlist = extlist,
                               pattern = pattern,
                               recursive = recursive,
                               pointer = TRUE,
                               cores = cores,
                               zap = zap)

  if(!is.null(keyvalues_out)){
    which_overlap = which(foreach(i = 1:length(image_list), .combine='c')%dopar%{
      Rwcs_overlap(image_list[[i]]$keyvalues, keyvalues_ref = keyvalues_out, buffer=0)
    })
    image_list = image_list[which_overlap]
  }else{
    which_overlap = 1:length(image_list)
    temp_overlap = 1:length(image_list)
  }

  if(!is.null(keyvalues_out)){
    if(isTRUE(keyvalues_out$ZIMAGE)){
      NAXIS1 = keyvalues_out$ZNAXIS1
      NAXIS2 = keyvalues_out$ZNAXIS2
    }else{
      NAXIS1 = keyvalues_out$NAXIS1
      NAXIS2 = keyvalues_out$NAXIS2
    }
  }else{
    NAXIS1 = dim(image_list[[1]])[1]
    NAXIS2 = dim(image_list[[1]])[2]
  }

  stack_grid = expand.grid(seq(1L,NAXIS1,by=chunk), seq(1L,NAXIS2,by=chunk))
  stack_grid[,3] = stack_grid[,1] + chunk
  stack_grid[stack_grid[,3] > NAXIS1,3] = NAXIS1
  stack_grid[,4] = stack_grid[,2] + chunk
  stack_grid[stack_grid[,4] > NAXIS2,4] = NAXIS2

  stack_med = foreach(i = 1:dim(stack_grid)[1])%dopar%{
  #for(i in 1:dim(stack_grid)[1]){
    message('Stacking sub region ',i,' of ',dim(stack_grid)[1])

    xsub = as.integer(stack_grid[i, c(1,3)])
    ysub = as.integer(stack_grid[i, c(2,4)])

    if(!is.null(keyvalues_out)){
      keyvalues_sub = Rwcs_keyvalues_sub(keyvalues_out, xsub=xsub, ysub=ysub)

      temp_overlap = which(foreach(j = 1:length(image_list), .combine = 'c')%dopar%{
        Rwcs_overlap(image_list[[j]]$keyvalues, keyvalues_sub, buffer=0)
      })
    }

    if(length(temp_overlap) == 0L){
      return(
        list(
          image = matrix(NA_real_, diff(range(xsub)) + 1L, diff(range(ysub)) + 1L),
          weight = matrix(0L, diff(range(xsub)) + 1L, diff(range(ysub)) + 1L)
        )
      )
    }

    image_list_cut = foreach(j = temp_overlap)%do%{ #not sure why this won't work in dopar... Rfits race conditions?
      if(!is.null(image_list[[j]]$keyvalues$XCUTLO) & useCUTLO){
        XCUTLO = image_list[[j]]$keyvalues$XCUTLO
      }else{
        XCUTLO = 1L #implictly assumed to be aligned images if XCUTLO is missing
      }

      if(!is.null(image_list[[j]]$keyvalues$YCUTLO) & useCUTLO){
        YCUTLO = image_list[[j]]$keyvalues$YCUTLO
      }else{
        YCUTLO = 1L #implictly assumed to be aligned images if YCUTLO is missing
      }

      xrange = c(1,(diff(range(xsub)) + 1L)) + (xsub[1] - XCUTLO)
      yrange = c(1,(diff(range(ysub)) + 1L)) + (ysub[1] - YCUTLO)
      return(imager::as.cimg(image_list[[j]][xrange, yrange]$imDat))
    }

    image = as.matrix(imager::parmed(image_list_cut, na.rm=TRUE))

    if(doweight){
      weight_list = foreach(j = 1:length(image_list_cut))%do%{
        return(imager::as.cimg(!is.na(image_list_cut[[j]])))
      }

      weight = as.matrix(imager::add(weight_list))
    }else{
      weight = NULL
    }

    return(
      list(image = image, weight = weight)
    )
  }

  image_out = matrix(0, NAXIS1, NAXIS2)
  if(doweight){
    weight_out = matrix(0L, NAXIS1, NAXIS2)
  }else{
    weight_out = NULL
  }

  for(i in 1:dim(stack_grid)[1]){
    image_out[stack_grid[i,1]:stack_grid[i,3], stack_grid[i,2]:stack_grid[i,4]] = stack_med[[i]]$image
    if(doweight){
      weight_out[stack_grid[i,1]:stack_grid[i,3], stack_grid[i,2]:stack_grid[i,4]] = stack_med[[i]]$weight
    }
  }

  if(!is.null(keyvalues_out)){
    keyvalues_out$EXTNAME = 'image'
    keyvalues_out$MAGZERO = image_list[[1]]$keyvalues$MAGZERO
    keyvalues_out$R_VER = R.version$version.string
    keyvalues_out$PANE_VER = as.character(packageVersion('ProPane'))
    keyvalues_out$RWCS_VER = as.character(packageVersion('Rwcs'))

    image_out = Rfits_create_image(image=image_out,
                                          keyvalues=keyvalues_out,
                                          keypass=FALSE,
                                          history='Stacked with Rwcs_stack')

    if(doweight){
      keyvalues_out$EXTNAME = 'weight'
      keyvalues_out$MAGZERO = NULL
      weight_out = Rfits_create_image(image=weight_out,
                                      keyvalues=keyvalues_out,
                                      keypass=FALSE)
    }
  }

  time_taken = proc.time()[3] - timestart
  message('Time taken: ',signif(time_taken,4),' seconds')

  output = list(
    image = image_out,
    weight = weight_out,
    which_overlap = which_overlap,
    time = time_taken,
    Nim = length(which_overlap)
  )

  class(output) = "ProPane"
  return(invisible(output))
}
