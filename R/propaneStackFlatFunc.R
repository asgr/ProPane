propaneStackFlatFunc = function(image_list=NULL, imager_func=NULL, na.rm=TRUE, weights=NULL,
                                increasing=TRUE, probs=0.5, ondisk=FALSE, cores = floor(detectCores()/2),
                                multitype = 'fork', chunk=1e3){

  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it from CRAN', call. = FALSE)
  }

  if(multitype=='fork'){
    registerDoParallel(cores=cores)
  }else if(multitype=='cluster'){
    registerDoParallel(cl=cores)
  }

  j = NULL

  for(i in 1:length(image_list)){
    if(ondisk){
      if(!inherits(image_list[[i]], 'Rfits_pointer')){
        stop('image_list[[',i,']] is not a Rfits_pointer')
      }

      if(i == 1L){
        NAXIS1 = dim(image_list[[i]])[1]
        NAXIS2 = dim(image_list[[i]])[2]
      }else{
        if(dim(image_list[[i]])[1] != NAXIS1){
          stop('image_list[[',i,']] NAXIS1 is not the same (',dim(image_list[[i]])[1],'versus ',NAXIS1,')')
        }
        if(dim(image_list[[i]])[2] != NAXIS2){
          stop('image_list[[',i,']] NAXIS2 is not the same (',dim(image_list[[i]])[2],'versus ',NAXIS2,')')
        }
      }
    }else{
      class(image_list) = 'list'
      if(inherits(image_list[[i]], 'Rfits_image')){
        image_list[[i]] = imager::as.cimg(image_list[[i]]$imDat)
      }else if(inherits(image_list[[i]], 'Rfits_pointer')){
        image_list[[i]] = imager::as.cimg(image_list[[i]][header=FALSE])
      }else if(!imager::is.cimg(image_list[[i]])){
        image_list[[i]] = imager::as.cimg(image_list[[i]])
      }
    }
  }

  if(ondisk){
    if(NAXIS1 %% chunk == 1L){
      stop('Chunk leaves single pixel sub-array, change chunk size!')
    }

    if(NAXIS2 %% chunk == 1L){
      stop('Chunk leaves single pixel sub-array, change chunk size!')
    }

    stack_grid = expand.grid(seq(1L,NAXIS1,by=chunk), seq(1L,NAXIS2,by=chunk))
    stack_grid[,3] = stack_grid[,1] + chunk
    stack_grid[stack_grid[,3] > NAXIS1,3] = NAXIS1
    stack_grid[,4] = stack_grid[,2] + chunk
    stack_grid[stack_grid[,4] > NAXIS2,4] = NAXIS2

    stack_sub = foreach(i = 1:dim(stack_grid)[1])%dopar%{
      message('Stacking sub region ',i,' of ',dim(stack_grid)[1])

      xsub = as.integer(stack_grid[i, c(1,3)])
      ysub = as.integer(stack_grid[i, c(2,4)])

      image_list_cut = foreach(j = 1:length(image_list))%do%{ #not sure why this won't work in dopar... Rfits race conditions?
        #header=FALSE since image_list[[j]] will be Rfits_pointer and we don't need header
        return(imager::as.cimg(image_list[[j]][xsub, ysub, header=FALSE]))
      }

      #Now we run our normal stacking just on the sub region
      output = propaneStackFlatFunc(
          image_list = image_list_cut,
          imager_func = imager_func,
          na.rm = na.rm,
          weights = weights,
          increasing = increasing,
          probs = probs,
          ondisk = FALSE,
      )

      return(output)
    }

    image_stack = matrix(0, NAXIS1, NAXIS2)
    weight_stack = matrix(0L, NAXIS1, NAXIS2)

    for(i in 1:dim(stack_grid)[1]){
      image_stack[stack_grid[i,1]:stack_grid[i,3], stack_grid[i,2]:stack_grid[i,4]] = stack_sub[[i]]$image
      weight_stack[stack_grid[i,1]:stack_grid[i,3], stack_grid[i,2]:stack_grid[i,4]] = stack_sub[[i]]$weight
    }

    return(invisible(list(image=image_stack, weight=weight_stack)))
  }else{
    if(na.rm == TRUE){
      weight_list = list()
      for(i in 1:length(image_list)){
        if(is.null(weights)){
          weight_list = c(weight_list, list(imager::as.cimg(!is.na(image_list[[i]]))))
        }else{
          weight_list = c(weight_list, list(imager::as.cimg(!is.na(image_list[[i]])*weights[i])))
        }
      }

      weight_stack = as.matrix(imager::add(weight_list))
      rm(weight_list)
    }else{
      weight_stack = NULL
    }

    if(is.null(imager_func)){
      imager_func = imager::average
    }

    if(is.function(imager_func)){
      if('w' %in% names(formals(imager_func))){
        if(is.null(weights)){
          stop('weights must be specified for imager wsum')
        }
        image_stack = as.matrix(imager_func(image_list, w=weights, na.rm=na.rm)) #only relevant for wsum
      }else if('prob' %in% names(formals(imager_func))){
        image_stack = as.matrix(imager_func(image_list, prob=probs, na.rm=na.rm))
      }else if('increasing' %in% names(formals(imager_func))){
        image_stack = imager_func(image_list, increasing=TRUE)
      }else{
        if('na.rm' %in% names(formals(imager_func))){
          image_stack = as.matrix(imager_func(image_list, na.rm=na.rm))
        }else{
          image_stack = as.matrix(imager_func(image_list))
        }
      }
    }else{
      #special cases:
      imager_func = tolower(imager_func)

      if(imager_func == 'quantile' | imager_func == 'quan' ){
        im_dim = dim(image_list[[1]])
        temp_mat = matrix(as.numeric(unlist(image_list)), nrow=length(image_list), byrow = TRUE)

        if(is.null(weights)){
          temp_weight = NULL
        }else{
          temp_weight = matrix(as.numeric(unlist(weights)), nrow=length(image_list_cut), byrow = TRUE)
        }

        temp_out = ProUtils::quan_wt_mat_col(temp_mat, probs=probs, wt=temp_weight)

        image_stack = list()

        if(length(probs) == 1){
          image_stack = matrix(temp_out, im_dim[1], im_dim[2])
        }else{
          for(k in 1:length(probs)){
            image_stack = c(image_stack, list(matrix(temp_out[,k], im_dim[1], im_dim[2])))
          }
        }

        #old code (slow!)
        # if(na.rm){
        #   image_stack = list()
        #   for(i in 1:length(prob)){
        #     temp_out = apply(temp_mat, MARGIN=2, FUN=quantile, probs=prob[i], na.rm=na.rm)
        #     image_stack = c(image_stack, list(matrix(temp_out, im_dim[1], im_dim[2])))
        #   }
        # }else{
        #   if(!requireNamespace("Rfast2", quietly = TRUE)){
        #     stop('The Rfast2 package is needed for this function to work. Please install it from CRAN', call. = FALSE)
        #   }
        #   temp_out = Rfast2::colQuantile(temp_mat, probs=prob, parallel=TRUE)
        #   image_stack = list()
        #   for(i in 1:length(prob)){
        #     image_stack = c(image_stack, list(matrix(temp_out[i,], im_dim[1], im_dim[2])))
        #   }
        # }
      }
    }

    return(invisible(list(image=image_stack, weight=weight_stack)))
  }
}

propaneStackFlatMed = function(image_list=NULL, na.rm=TRUE, ondisk=FALSE,
                               cores = floor(detectCores()/2), multitype = 'fork', chunk=1e3){

  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for stacking to work. Please install from GitHub asgr/Rfits.', call. = FALSE)
  }

  output = propaneStackFlatFunc(
    image_list = image_list,
    imager_func = imager::parmed,
    na.rm = na.rm,
    ondisk = ondisk,
    cores = cores,
    multitype = multitype,
    chunk = chunk
  )

  return(output)

}
