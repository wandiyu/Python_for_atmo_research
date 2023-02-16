def GRID_3d(lon, lat, vert, lonGrid, latGrid, vertGrid, varin):
    
    # grid 3-d data: lon, lat, vertical coordinate 
    # each observation has different lat,lon, and pressure/altitude, example: SABER 
    # first bin observations into lat-lon grids
    # then linear interpolate all observations in bin into vertical grid 
    
    import numpy as np
    from scipy.interpolate import interp1d
    lon = lon.flatten()
    lat = lat.flatten()
    vert = vert.flatten()
    varin = [v.flatten() for v in varin]
    
    nlatGrid  = len(latGrid)                           # the number of latitude grids
    nlonGrid  = len(lonGrid)                           # the number of longitude grids
    nvertGrid = len(vertGrid)                          # the number of altitude grids
    nvar = len(varin)                                  # the number of variables
    
    latStep   = latGrid[1] - latGrid[0]                # calculate the latstep for the regular latitudes
    lonStep   = lonGrid[1] - lonGrid[0]                # calculate the lonstep for the regular longitude

    varout = np.ones((nvar,nvertGrid,nlatGrid,nlonGrid))*np.nan
    #grid the variable
    for j in range(nlatGrid):
        for i in range(nlonGrid):
            lon_low = lonGrid[i]-lonStep/2.0 ; lon_high = lonGrid[i]+lonStep/2.0
            lat_low = latGrid[j]-latStep/2.0 ; lat_high = latGrid[j]+latStep/2.0
            iH = np.where(((lon >= lon_low) & (lon < lon_high)) & ((lat >= lat_low) & (lat < lat_high)))[0]
            if len(iH) != 0:
                for i_var in range(nvar):
                    x = varin[i_var][iH]
                    if len(x[np.isfinite(x)])>=2:
                        finterp = interp1d(vert[iH],varin[i_var][iH],fill_value=np.nan,
                                bounds_error=False)
                        varout[i_var,:,j,i] = finterp(vertGrid)
    return varout



def GRID_2d(lon, lat, lonGrid, latGrid, varin):
    
    # grid 2-d data: lon, lat
    # each observation has different lat,lon, but same pressure/altitude, example: MLS 
    # bin observations into lat-lon grids and then average
    
    import numpy as np
    
    nlatGrid  = len(latGrid)                           # the number of latitude grids
    nlonGrid  = len(lonGrid)                           # the number of longitude grids
    nvar = len(varin)                                  # the number of variables
    
    latStep   = latGrid[1] - latGrid[0]                # calculate the latstep for the regular latitudes
    lonStep   = lonGrid[1] - lonGrid[0]                # calculate the lonstep for the regular longitude
    var_shape = varin[0].shape
    varout = np.zeros((nvar,var_shape[1],nlatGrid,nlonGrid))
    #grid the variable
    for j in range(nlatGrid):
        for i in range(nlonGrid):
            lon_low = lonGrid[i]-lonStep/2.0 ; lon_high = lonGrid[i]+lonStep/2.0
            lat_low = latGrid[j]-latStep/2.0 ; lat_high = latGrid[j]+latStep/2.0
            iH = np.where(((lon >= lon_low) & (lon < lon_high)) & ((lat >= lat_low) & (lat < lat_high)))
            if len(iH[0]) != 0:
                for k in range(var_shape[1]):
                    for i_var in range(nvar):
                        varout[i_var,k,j,i] = np.nanmean(varin[i_var][iH[0],k])
    return varout


def GRID_3d_sameloc(lon, lat, vert, lonGrid, latGrid, vertGrid, varin):
    
    # grid 3-d data: lon, lat, vertical coordinate
    # each observation has the same lat,lon, but different pressure/altitude, example: HALOE  
    # linear interpolate into vertical coordinates  
    
    
    import numpy as np
    from scipy.interpolate import interp1d
    import itertools
    
    nlatGrid  = len(latGrid)                           # the number of latitude grids
    nlonGrid  = len(lonGrid)                           # the number of longitude grids
    nvertGrid = len(vertGrid)                          # the number of altitude grids
    nvar = len(varin)                                  # the number of variables
    
    latStep   = latGrid[1] - latGrid[0]                # calculate the latstep for the regular latitudes
    lonStep   = lonGrid[1] - lonGrid[0]                # calculate the lonstep for the regular longitude

    varout = np.zeros((nvar,nvertGrid,nlatGrid,nlonGrid))
    #grid the variable
    for j in range(nlatGrid):
        for i in range(nlonGrid):
            lon_low = lonGrid[i]-lonStep/2.0 ; lon_high = lonGrid[i]+lonStep/2.0
            lat_low = latGrid[j]-latStep/2.0 ; lat_high = latGrid[j]+latStep/2.0
            iH = np.where(((lon >= lon_low) & (lon < lon_high)) & ((lat >= lat_low) & (lat < lat_high)))[0]
            if len(iH) != 0:
                for i_var in range(nvar):
                    y = [varin[i_var][iHH] for iHH in iH]
                    y = np.array(list(itertools.chain.from_iterable(y)))
                    x = [vert[iHH] for iHH in iH]
                    x = np.array(list(itertools.chain.from_iterable(x)))
                    if len(y[np.isfinite(y)])>=2:
                        finterp = interp1d(x,y,fill_value=np.nan,
                                bounds_error=False)
                        varout[i_var,:,j,i] = finterp(vertGrid)
    return varout


    
