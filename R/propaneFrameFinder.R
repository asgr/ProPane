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

  # MatchInfo_Cor_BL = coordmatchsing(RAcen,
  #                                Deccen,
  #                                FrameInfo[,c('corner_BL_RA','corner_BL_Dec')], rad=rad, radunit=radunit)
  #
  # MatchInfo_Cor_TL = coordmatchsing(RAcen,
  #                                   Deccen,
  #                                   FrameInfo[,c('corner_TL_RA','corner_TL_Dec')], rad=rad, radunit=radunit)
  #
  # MatchInfo_Cor_TR = coordmatchsing(RAcen,
  #                                   Deccen,
  #                                   FrameInfo[,c('corner_TR_RA','corner_TR_Dec')], rad=rad, radunit=radunit)
  #
  # MatchInfo_Cor_BR = coordmatchsing(RAcen,
  #                                   Deccen,
  #                                   FrameInfo[,c('corner_BR_RA','corner_BR_Dec')], rad=rad, radunit=radunit)

  # Pointers = Rfits_make_list(filelist = filelist,
  #                            dirlist = dirlist,
  #                            extlist = extlist,
  #                            cores = cores)

  # registerDoParallel(cores=cores)

  # point = NULL

  # In_Frame = FrameInfo[,'min_RA'] < RAcen & FrameInfo[,'max_RA'] > RAcen &
  #            FrameInfo[,'min_Dec'] < Deccen & FrameInfo[,'max_Dec'] > Deccen

  # GoodID = c(MatchInfo_Cen$ID,
  #            MatchInfo_Cor_BL$ID,
  #            MatchInfo_Cor_TL$ID,
  #            MatchInfo_Cor_TR$ID,
  #            MatchInfo_Cor_BR$ID,
  #            which(In_Frame)
  #            )

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
  #FrameInfo = cbind(FrameInfo, sep = MatchInfo$sep)

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
  legend('topleft', legend = paste0('Nmatch = ', Nmatch))

  return(invisible(FrameInfo))
}
