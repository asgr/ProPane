propaneDistortPixscale = function(keyvalues, useraw=TRUE, unit='asec', ...){
  im_dim = dim(keyvalues)
  pixgrid = expand.grid(1:im_dim[1],1:im_dim[2])

  if(useraw){
    header = Rfits_keyvalues_to_raw(keyvalues)
  }else{
    header = NULL
  }

  BL = Rwcs_p2s(pixgrid[,1] - 0.5, pixgrid[,2] - 0.5, keyvalues=keyvalues, header=header, pixcen='FITS', ...)
  BR = Rwcs_p2s(pixgrid[,1] + 0.5, pixgrid[,2] - 0.5, keyvalues=keyvalues, header=header, pixcen='FITS', ...)
  TL = Rwcs_p2s(pixgrid[,1] - 0.5, pixgrid[,2] + 0.5, keyvalues=keyvalues, header=header, pixcen='FITS', ...)

  BL[,1] = BL[,1] * cos(BL[,2]*pi/180)
  BR[,1] = BR[,1] * cos(BR[,2]*pi/180)
  TL[,1] = TL[,1] * cos(TL[,2]*pi/180)

  scale = 0.7071068*sqrt((BR[,1] - BL[,1])^2 + (BR[,2] - BL[,2])^2 + (TL[,1] - BL[,1])^2 + (TL[,2] - BL[,2])^2)

  if(unit=='deg'){
    #do nothing
  }else if(unit == 'asec'){
    scale = scale*3600
  }else if(unit == 'amin'){
    scale = scale*60
  }else if(unit=='rad'){
    scale = scale*(pi/180)
  }else{
    message('Not a valid unit, must be one of asec2 / amin2 / deg2 / rad2 / str')
  }

  output = Rfits_create_image(matrix(scale, im_dim[1], im_dim[2]), keyvalues=keyvalues)

  return(output)
}

propaneDistortPixarea = function(keyvalues, useraw=TRUE, unit='asec2', ...){
  im_dim = dim(keyvalues)
  pixgrid = expand.grid(1:im_dim[1],1:im_dim[2])

  if(useraw){
    header = Rfits_keyvalues_to_raw(keyvalues)
  }else{
    header = NULL
  }

  BL = Rwcs_p2s(pixgrid[,1] - 0.5, pixgrid[,2] - 0.5, keyvalues=keyvalues, header=header, pixcen='FITS', ...)
  BR = Rwcs_p2s(pixgrid[,1] + 0.5, pixgrid[,2] - 0.5, keyvalues=keyvalues, header=header, pixcen='FITS', ...)
  TL = Rwcs_p2s(pixgrid[,1] - 0.5, pixgrid[,2] + 0.5, keyvalues=keyvalues, header=header, pixcen='FITS', ...)

  BL[,1] = BL[,1] * cos(BL[,2]*pi/180)
  BR[,1] = BR[,1] * cos(BR[,2]*pi/180)
  TL[,1] = TL[,1] * cos(TL[,2]*pi/180)

  area = sqrt((BR[,1] - BL[,1])^2 + (BR[,2] - BL[,2])^2)*sqrt((TL[,1] - BL[,1])^2 + (TL[,2] - BL[,2])^2) #in deg

  if(unit=='deg2'){
    #do nothing
  }else if(unit == 'asec2'){
    area = area*3600^2
  }else if(unit == 'amin2'){
    area = area*60^2
  }else if(unit=='rad2' | unit=='str'){
    area = area*(pi/180)^2
  }else{
    message('Not a valid unit, must be one of asec2 / amin2 / deg2 / rad2 / str')
  }

  output = Rfits_create_image(matrix(area, im_dim[1], im_dim[2]), keyvalues=keyvalues)

  return(output)
}
