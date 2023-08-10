propaneStackFlatFunc = function(image_list=NULL, imager_func=NULL, na.rm=TRUE, weights=NULL,
                                increasing=TRUE, probs=c(0.159, 0.5, 0.841), parallel=TRUE){

  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it from CRAN', call. = FALSE)
  }

  for(i in 1:length(image_list)){
    if(inherits(image_list[[i]], 'Rfits_image')){
      image_list[[i]] = imager::as.cimg(image_list[[i]]$imDat)
    }
    if(inherits(image_list[[i]], 'Rfits_pointer')){
      image_list[[i]] = imager::as.cimg(image_list[[i]][header=FALSE])
    }
    if(!imager::is.cimg(image_list[[i]])){
      image_list[[i]] = imager::as.cimg(image_list[[i]])
    }
  }

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
    if('na.rm' %in% names(formals(imager_func))){
      if('w' %in% names(formals(imager_func))){
        image_stack = as.matrix(imager_func(image_list, w=weights, na.rm=na.rm))
      }else{
        image_stack = as.matrix(imager_func(image_list, na.rm=na.rm))
      }
    }else{
      if('increasing' %in% names(formals(imager_func))){
        image_stack = imager_func(image_list, increasing=increasing)
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

      if(na.rm){
        image_stack = list()
        for(i in 1:length(probs)){
          temp_out = apply(temp_mat, MARGIN=2, FUN=quantile, probs=probs[i], na.rm=TRUE)
          image_stack = c(image_stack, list(matrix(temp_out, im_dim[1], im_dim[2])))
        }
      }else{
        if(!requireNamespace("Rfast2", quietly = TRUE)){
          stop('The Rfast2 package is needed for this function to work. Please install it from CRAN', call. = FALSE)
        }
        temp_out = Rfast2::colQuantile(temp_mat, probs=probs, parallel=parallel)
        image_stack = list()
        for(i in 1:length(probs)){
          image_stack = c(image_stack, list(matrix(temp_out[i,], im_dim[1], im_dim[2])))
        }
      }
    }
  }

  return(invisible(list(image=image_stack, weight=weight_stack)))
}
