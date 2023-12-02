propaneFrameFinder = function(filelist = NULL,
                            dirlist = NULL,
                            extlist = 1,
                            RAcen = 180,
                            Deccen = 0,
                            rad = 1,
                            cores = 1,
                            plot = TRUE,
                            proj = TRUE,
                            border = 'red',
                            col = hsv(alpha=0.2),
                            bg = par("bg"),
                            ...){

  FrameInfo = Rfits_key_scan(filelist = filelist,
                             dirlist = dirlist,
                             extlist = extlist,
                             cores = cores,
                             get_dim = TRUE,
                             get_centre = TRUE,
                             get_corners = TRUE,
                             get_pixscale = TRUE,
                             data.table = FALSE,
                             ...)

  MatchInfo_Cen = coordmatchsing(RAcen,
                           Deccen,
                           FrameInfo[,c('centre_RA','centre_Dec')], rad=rad + FrameInfo$pixscale*sqrt(FrameInfo$dim_1^2 + FrameInfo$dim_2^2)/2/3600, radunit='deg')


  GoodID = MatchInfo_Cen$ID

  GoodID = GoodID[!is.na(GoodID)]

  if(length(GoodID) == 0L){
    message('No files overlap with target region!')
    return(NULL)
  }else{
    GoodID = sort(unique(GoodID))

    Nmatch = length(GoodID)
  }

  FrameInfo = FrameInfo[GoodID,]

  longlim = RAcen + c(-rad,rad)/cos(Deccen*pi/180)
  latlim = Deccen + c(-rad,rad)

  if(plot == FALSE){
    return(invisible(FrameInfo))
  }

  if(proj){
    magproj(RAcen, Deccen,
            type = 'p',
            longlim = longlim,
            latlim = latlim,
            centre = c(RAcen, Deccen),
            labloc = c(longlim[2] + rad/cos(Deccen*pi/180)/20, latlim[1] -  rad/20),
            pch = 4, cex = 5,
            fliplong = TRUE
    )
    mtext('RA / deg', side = 1)
    mtext('Dec / deg', side = 2)

    for(i in 1:dim(FrameInfo)[1]){
      magproj(as.numeric(FrameInfo[i,c('corner_BL_RA','corner_TL_RA','corner_TR_RA','corner_BR_RA')]),
              as.numeric(FrameInfo[i,c('corner_BL_Dec','corner_TL_Dec','corner_TR_Dec','corner_BR_Dec')]),
              type = 'pl',
              border = border,
              col = col,
              add = TRUE
      )
    }
  }else{
    magplot(RAcen, Deccen,
            xlim = rev(longlim),
            ylim = latlim,
            asp = 1/cos(Deccen*pi/180),
            pch = 4, cex = 5,
            xlab = 'RA / deg',
            ylab = 'Dec / deg'
            )

    for(i in 1:dim(FrameInfo)[1]){
      polygon(as.numeric(FrameInfo[i,c('corner_BL_RA','corner_TL_RA','corner_TR_RA','corner_BR_RA')]),
              as.numeric(FrameInfo[i,c('corner_BL_Dec','corner_TL_Dec','corner_TR_Dec','corner_BR_Dec')]),
              border = border,
              col = col
              )
    }
  }
  legend('topleft', legend = paste('Nmatch:', Nmatch), bg=bg)

  return(invisible(FrameInfo))
}

propaneGenWCS = function(filelist = NULL, dirlist = NULL, image_list = NULL, rotation = 'North',
                         CRVAL1 = NULL, CRVAL2 = NULL, pixscale = NULL,
                         NAXIS1 = NULL, NAXIS2 = NULL, CRPIX1 = NULL, CRPIX2 = NULL,
                         CTYPE1 = "RA---TAN", CTYPE2 = "DEC--TAN",
                         CD1_1 = NULL, CD1_2 = NULL, CD2_1 = NULL, CD2_2 = NULL,
                         CUNIT1 = "deg", CUNIT2 = "deg", ...){
  
    info = Rfits_key_scan(filelist = filelist,
                          dirlist = dirlist,
                          image_list = image_list,
                          get_centre = TRUE,
                          get_rotation = TRUE,
                          get_corners = TRUE,
                          get_extremes = TRUE,
                          get_pixscale = TRUE,
                          ...)

  if(is.null(CRVAL1)){
    CRVAL1 = (min(info$min_RA) + max(info$max_RA))/2
  }
  if(is.null(CRVAL2)){
    CRVAL2 = (min(info$min_Dec) + max(info$max_Dec))/2
  }
  if(is.null(pixscale)){
    pixscale = min(info$pixscale)
  }
  if(is.null(NAXIS1)){
    NAXIS1 = ceiling(3600*diff(range(info$min_RA, info$max_RA)) * cos(max(abs(info$min_Dec), abs(info$max_Dec))*pi/180)/pixscale)
  }
  if(is.null(NAXIS2)){
    NAXIS2 = ceiling(3600*diff(range(info$min_Dec, info$max_Dec))/pixscale)
  }
  if(is.null(CRPIX1)){
    CRPIX1 = NAXIS1/2 + 0.5
  }
  if(is.null(CRPIX2)){
    CRPIX2 = NAXIS2/2 + 0.5
  }
  if(is.null(CD1_1)){
    CD1_1 = -pixscale/3600
  }
  if(is.null(CD1_2)){
    CD1_2 = 0
  }
  if(is.null(CD2_1)){
    CD2_1 = 0
  }
  if(is.null(CD2_2)){
    CD2_2 = pixscale/3600
  }

  keyvalues = Rwcs_setkeyvalues(
    CRVAL1 = CRVAL1,
    CRVAL2 = CRVAL2,
    NAXIS1 = NAXIS1,
    NAXIS2 = NAXIS2,
    CRPIX1 = CRPIX1,
    CRPIX2 = CRPIX2,
    CTYPE1 = CTYPE1,
    CTYPE2 = CTYPE2,
    CUNIT1 = CUNIT1,
    CUNIT2 = CUNIT2,
    CD1_1 = CD1_1,
    CD1_2 = CD1_2,
    CD2_1 = CD2_1,
    CD2_2 = CD2_2
  )
  
  if(is.character(rotation)){
    rotation = tolower(rotation)
    if(rotation == 'north'){rotation = 0}
    if(rotation == 'east'){rotation = 90}
    if(rotation == 'south'){rotation = 180}
    if(rotation == 'west'){rotation = 270}
    
    if(rotation == 'get'){
      temp_rot = info$rotation_North
      temp_rot = temp_rot %% 90
      temp_rot[temp_rot > 45] = temp_rot[temp_rot > 45] - 90
      rotation = mean(temp_rot, na.rm=TRUE)
    }
  }
  
  if(rotation != 0){
    keyvalues = propaneWCSmod(keyvalues, delta_rot = -rotation)
  }
  
  RAcorners = c(info$corner_BL_RA, info$corner_BR_RA, info$corner_TL_RA, info$corner_TR_RA)
  Deccorners = c(info$corner_BL_Dec, info$corner_BR_Dec, info$corner_TL_Dec, info$corner_TR_Dec)
  
  xycorners = Rwcs_s2p(RAcorners, Deccorners, keyvalues=keyvalues)
  
  xlo = ceiling(min(xycorners[,'x'], na.rm=TRUE))
  xhi = ceiling(max(xycorners[,'x'], na.rm=TRUE))
  ylo = ceiling(min(xycorners[,'y'], na.rm=TRUE))
  yhi = ceiling(max(xycorners[,'y'], na.rm=TRUE))
  
  if(xlo != 1L){
    keyvalues$CRPIX1 = keyvalues$CRPIX1 - xlo + 1L
  }
  if(ylo != 1L){
    keyvalues$CRPIX2 = keyvalues$CRPIX2 - ylo + 1L
  }
  
  keyvalues$NAXIS1 = xhi - xlo + 1L
  keyvalues$NAXIS2 = yhi - ylo + 1L

  return(keyvalues)
}
