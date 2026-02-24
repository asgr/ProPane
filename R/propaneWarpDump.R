propaneWarpDump = function(image_list=NULL, magzero_in=0, magzero_out=23.9, keyvalues_out=NULL,
                           dim_out=NULL, cores=floor(detectCores()/2), cores_warp=1,  multitype='fork',
                           keepcrop=TRUE, dump_dir=tempdir(), stub='image', ...){

  if(!requireNamespace("Rfits", quietly = TRUE)){
    stop('The Rfits package is needed for stacking to work. Please install from GitHub asgr/Rfits.', call. = FALSE)
  }

  if(!requireNamespace("Rwcs", quietly = TRUE)){
    stop('The Rfits package is needed for stacking to work. Please install from GitHub asgr/Rfits.', call. = FALSE)
  }

  dump_dir = path.expand(dump_dir)
  if(dir.exists(dump_dir) == FALSE){
    dir.create(dump_dir, recursive = TRUE)
  }
  message('Frames being dumped to ', dump_dir)

  if(multitype=='fork'){
    registerDoParallel(cores=cores)
  }else if(multitype=='cluster'){
    registerDoParallel(cl=cores)
  }

  Nim = length(image_list)

  if(length(magzero_in) == 1){
    magzero_in = rep(magzero_in, Nim)
  }

  if(is.null(keyvalues_out)){
    keyvalues_out = image_list[[1]]$keyvalues
  }

  if(is.null(keyvalues_out)){
    stop('Need keyvalues out!')
  }

  if(isTRUE(keyvalues_out$ZIMAGE)){
    dim_im = c(keyvalues_out$ZNAXIS1, keyvalues_out$ZNAXIS2)
  }else{
    dim_im = c(keyvalues_out$NAXIS1, keyvalues_out$NAXIS2)
  }

  # Check all supplied frames are in WCS:

  i = NULL
  which_overlap = which(foreach(i = 1:Nim, .combine='c')%dopar%{
    Rwcs_overlap(image_list[[i]]$keyvalues, keyvalues_ref = keyvalues_out, buffer=0)
  })

  Ncheck = length(which_overlap)

  if(Ncheck == 0){
    message('No frames exist within target WCS!')
    return(NULL)
  }

  if(Nim > Ncheck){
    message('Only ', Ncheck, ' of ', Nim, ' input frames overlap with target WCS!')

    Nim = Ncheck

    # Trim down all the inputs to just those overlapping with the target WCS:

    image_list = image_list[which_overlap]


    magzero_in = magzero_in[which_overlap]
  }

  zero_point_scale = 10^(-0.4*(magzero_in - magzero_out))

  message('Projecting Images 1 to ', Nim)

  foreach(i = 1:Nim)%dopar%{
    suppressMessages({
      temp_warp = propaneWarp(
        image_in = image_list[[i]],
        keyvalues_out = keyvalues_out,
        dim_out = dim_out,
        doscale = TRUE,
        dotightcrop = TRUE,
        keepcrop = keepcrop,
        warpfield_return = TRUE,
        cores = cores_warp,
        ...
      )

      if(zero_point_scale[i] != 1){
        temp_warp$imDat = temp_warp$imDat*zero_point_scale[i]
      }

      temp_warp$keyvalues$MAGZERO = magzero_out
      temp_warp$keycomments$MAGZERO = ""
      temp_warp$keynames['MAGZERO'] = "MAGZERO"
      Rfits_write_image(temp_warp, paste0(dump_dir,'/',stub,'_warp_',i,'.fits'))
    })
    return(NULL)
  }

  return(dump_dir)
}
