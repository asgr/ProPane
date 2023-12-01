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

propaneGenWCS = function(image_list = NULL, filelist = NULL, dirlist = NULL, CRVAL1 = NULL,
                          CRVAL2 = NULL, pixscale = NULL, NAXIS1 = NULL, NAXIS2 = NULL,
                          CRPIX1 = NULL, CRPIX2 = NULL, CTYPE1 = "RA---TAN", CTYPE2 = "DEC--TAN",
                          CUNIT1 = "deg", CUNIT2 = "deg", ...){
  if(!is.null(image_list)){
    image = NULL
    info = foreach(image = image_list, .combine='rbind')%do%{
      current_info = list()

      temp_cen = centre(image, ...)
      if(is.na(temp_cen[1])){
        temp_cen = rep(NA, 2)
      }
      current_info = c(current_info,
                       centre_RA = temp_cen[1], centre_Dec = temp_cen[2]
      )

      temp_cor = corners(image, ...)
      if(is.na(temp_cor[1])){
        temp_cor = matrix(NA, 4, 2)
      }
      current_info = c(current_info,
                       corner_BL_RA = temp_cor[1,1], corner_BL_Dec = temp_cor[1,2],
                       corner_TL_RA = temp_cor[2,1], corner_TL_Dec = temp_cor[2,2],
                       corner_TR_RA = temp_cor[3,1], corner_TR_Dec = temp_cor[3,2],
                       corner_BR_RA = temp_cor[4,1], corner_BR_Dec = temp_cor[4,2]
      )

      temp_ext = extremes(image, ...)
      if(is.na(temp_ext[1])){
        temp_ext = matrix(NA, 3, 2)
      }
      current_info = c(current_info,
                       min_RA = temp_ext[1,1], min_Dec = temp_ext[1,2],
                       max_RA = temp_ext[2,1], max_Dec = temp_ext[2,2],
                       range_RA = temp_ext[3,1], range_Dec = temp_ext[3,2]
      )

      temp_pixscale = pixscale(image, ...)
      current_info = c(current_info, pixscale = temp_pixscale)

      return(as.data.frame(current_info))
    }
  }else{
    info = Rfits_key_scan(filelist=filelist, dirlist=dirlist,
                                get_centre = TRUE,
                                get_corners = TRUE,
                                get_extremes = TRUE,
                                get_pixscale = TRUE,
                                ...)
  }

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
    NAXIS1/2 + 0.5
  }
  if(is.null(CRPIX2)){
    NAXIS2/2 + 0.5
  }

  keyvalues = Rwcs_setkeyvalues(
    CRVAL1 = CRVAL1,
    CRVAL2 = CRVAL2,
    pixscale = pixscale,
    NAXIS1 = NAXIS1,
    NAXIS2 = NAXIS2,
    CRPIX1 = CRPIX1,
    CRPIX2 = CRPIX2,
    CTYPE1 = CTYPE1,
    CTYPE2 = CTYPE2,
    CUNIT1 = CUNIT1,
    CUNIT2 = CUNIT2
  )

  return(keyvalues)
}
