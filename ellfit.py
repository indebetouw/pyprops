def ellfit(x,y,wt=None):
    import pylab as pl
    # Calculate the best fit ellipse for an X and Y distribution, allowing
    # for weighting.
    # OUTPUTS:
    #   MAJOR - major axis in same units as x and y
    #   MINOR - minor axis in same units as x and y
    #   POSANG - the position angle CCW from the X=0 line of the coordinates 
    #
    #   Adam: The intensity weighted major and minor values are equal to the
    #   second moment. For equal weighting by pixel (of the sort that
    #   might be done for blob analysis) the ellipse fit to the
    #   half-maximum area will have semimajor axis equal to 1./1.69536 the
    #   second moment. For the quarter maximum surface this is 1./1.19755.
    #
    #   but adam did not have the factor of 4 to turn eigenval into major axis 
    # 
    #   translation: if we run this with intensity weight, we get
    #   the second moment back (a sigma).  for flat weights i think he means
    #   the halfmax contour semimajor axis 
    
    if type(wt)==type(None):
        wt = x*0.0 + 1.0
        
    tot_wt = wt.sum()

# WEIGHTED X AND Y CENTERS
    x_ctr = (wt*x).sum()/tot_wt
    y_ctr = (wt*y).sum()/tot_wt

# BUILD THE MATRIX
    i11 = (wt*(x-x_ctr)**2).sum()/tot_wt
    i22 = (wt*(y-y_ctr)**2).sum()/tot_wt
    i12 = (wt*(x-x_ctr)*(y-y_ctr)).sum()/tot_wt
    mat = [[i11, i12], [i12, i22]]

# CATCH THE CASE OF ZERO DETERMINANT  
    if pl.det(mat)==0:
        return pl.nan,pl.nan,pl.nan

    if pl.any(pl.isnan(mat)):
        return pl.nan,pl.nan,pl.nan

# WORK OUT THE EIGENVALUES
    evals,evec = pl.eig(mat)

# PICK THE MAJOR AXIS
    absvals=pl.absolute(evals)
    major = absvals.max()
    maj_ind = pl.where(absvals==major)[0][0]
    major_vec = evec[maj_ind]
    min_ind = 1-maj_ind

# WORK OUT THE ORIENTATION OF THE MAJOR AXIS
    posang = pl.arctan2(major_vec[1], major_vec[0])

    # compared to the original idl code, this code is returning
    # pi-the desired angle, so:
#    posang=pl.pi-posang

#    if posang<0: posang = posang+pl.pi

# MAJOR AND MINOR AXIS SIZES
# turn into real half-max major/minor axis
    major = pl.sqrt(evals[maj_ind])*4.
    minor = pl.sqrt(evals[min_ind])*4.

    return major,minor,posang

#===========================================================================
