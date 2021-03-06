def extrap( xin, yin, targett = 0., 
            fast = True, square = False, weight = True):

# NAME:
#  EXTRAP
# PURPOSE:
#  To extrapolate a function T, Y(T) to a set value of T
#
# CALLING SEQUENCE:
#   value = EXTRAP(t, y, [targett = targett, /fast, /square, scatter =
#   scatter, /weight])
#
# INPUTS:
#   T, Y -- The indpendent and dependent variables for extrapolation. 
#   TARGETT -- The target T value for extrapolation
#
# KEYWORD PARAMETERS:
#   /FAST -- Subsample the data to speed up performance
#   /SQUARE -- Do a quadratic extrapolation
#   SCATTER -- The M.A.D. of the extrapolated moment values.
#   /WEIGHT -- Weight the data according to contour level.
#   COEFFS -- Set to a named variable to received the coefficients of
#             the fit.
#
# OUTPUTS:
#   VALUE -- The extrapolated value.
#
# MODIFICATION HISTORY:
#
#       Documented -- Fri Sep 2 16:47:24 2005, Erik Rosolowsky
#                     <erosolow@asgard.cfa.harvard.edu>
#

# scatter = !values.f_nan

    import matplotlib.pyplot as pl
    import numpy as np
    
# THE TARGET FOR THE EXTRAPOLATION
    if (targett>xin.min()) and (targett<xin.max()):
        d = np.absolute(xin - targett)
        z = np.where(d==d.min())[0]
        #if len(z)>1: # e.g. it goes down into noise, with targett=0
        return np.nanmedian(yin[z])

# CAN'T EXTRAPOLATE WITH LESS THAN 5 DATA POINTS
    if len(xin) < 5:
        return np.nan

# REVERSE THE (SORTED increasing) ARRAYS BEFORE EXTRAPOLATING
    xuse = xin[::-1]
    yuse = yin[::-1]

    import pylab as pl
# CULL OUT BAD VALUES
    goodind = np.where(np.isnan(yuse)==False)[0]
    if len(goodind)>0:
        yuse = yuse[goodind]
        xuse = xuse[goodind]
  
# TO DO THINGS FAST, WE DECIMATE THE ARRAY DOWN TO 250 ELEMENTS
    if fast and (len(xuse) > 250):
        useindex = np.int32(np.arange(250)*len(xuse)/250.)
        xuse = xuse[useindex]
        yuse = yuse[useindex]
  
# KEEP FITTING, INCLUDING PROGRESSIVELY MORE AND MORE OF THE CLOUD
    if not weight:

        coeffs = np.ones([2,len(xuse)])*np.nan
        if square:
            coeffs = np.ones([3,len(xuse)])*np.nan
    
        for i in np.arange(len(xuse)-3)+3:
            xfit = xuse[0:i]
            yfit = yuse[0:i]    
            n = len(xfit)
            
            xfits  = xfit.sum()
            yfits  = yfit.sum()
            xfit2s = (xfit**2).sum()
            yfit2s = (yfit**2).sum()
            
            if not square:
                M = [[n, xfits], 
                     [xfits, xfit2s]]
                covar = np.linalg.inv(M)
                coeffs[:, i] = np.dot(covar,[yfits, (xfit*yfit).sum()])
            else:
                xfit3s = (xfit**3).sum()
                M = [[n, xfits, xfit2s], 
                     [xfits, xfit2s, xfit3s], 
                     [xfit2s, xfit3s, (xfit**4).sum()]]
                covar = np.linalg.inv(M)
                coeffs[:, i] = np.dot(covar,[yfits, (xfit*yfit).sum(), (xfit**2*yfit).sum()])

        #   DO THE EXTRAPOLATION
        extrap_value = np.nanmedian(coeffs[0, :] + coeffs[1, :]*targett)
        def mad(data, axis=None):
            return np.nanmedian(np.absolute(data - np.nanmedian(data, axis)), axis)
        scatter =  mad(coeffs[0, :] + coeffs[1, :]*targett)
        if square:
            extrap_value = np.nanmedian(coeffs[0, :] + coeffs[1, :]*targett + 
                                     coeffs[2, :]*targett**2)  
            scatter = mad(coeffs[0, :] + coeffs[1, :]*targett + 
                          coeffs[2, :]*targett**2)  

    else:
        coeffs = np.ones(2)*np.nan
        if square:
            coeffs = np.ones(3)*np.nan

        xfit = xuse
        yfit = yuse
    
        # WEIGHT BY CONTOUR NUMBER (BIG = MORE WEIGHT)
        wfit = (np.arange(len(yfit))+1)[::-1]

        crosssum=(xfit*wfit).sum()
        if not square:            
            M = [[wfit.sum(), crosssum], 
                 [crosssum, (xfit**2*wfit).sum()]]
            covar = np.linalg.inv(M)
            coeffs = np.dot(covar,[(yfit*wfit).sum(), (xfit*yfit*wfit).sum()])
        else:
            crosssum2=(xfit**2*wfit).sum()
            crosssum3=(xfit**3*wfit).sum()
            M = [[wfit.sum(), crosssum, crosssum2], 
                 [crosssum, crosssum2, crosssum3], 
                 [crosssum2, crosssum3, (xfit**4*wfit).sum()]]
            covar = np.linalg.inv(M)
            coeffs = np.dot(covar,[(yfit*wfit).sum(), (xfit*yfit*wfit).sum(),
                                   (xfit**2*yfit*wfit).sum()])

        extrap_value = coeffs[0] + coeffs[1]*targett
        if square:
            extrap_value = extrap_value + coeffs[2]*targett**2
  
    return extrap_value
  
