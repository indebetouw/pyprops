def deconvolve_gauss(meas_maj,       # measured major axis
                     meas_min,       # measured minor axis
                     meas_pa,        # measured position angle
                     beam_maj,       # beam major axis
                     beam_min,       # beam minor axis
                     beam_pa):       # beam position angle

    #  ADAPTED FROM gaupar.for in MIRIAD via K. Sandstrom
    #  IO in radians
    meas_theta=meas_pa
    beam_theta=beam_pa

    import pylab as pl
    # MATH (FROM MIRIAD VIA K. SANDSTROM)
    alpha = \
        (meas_maj*pl.cos(meas_theta))**2 + (meas_min*pl.sin(meas_theta))**2 - \
        (beam_maj*pl.cos(beam_theta))**2 - (beam_min*pl.sin(beam_theta))**2

    beta = \
        (meas_maj*pl.sin(meas_theta))**2 + (meas_min*pl.cos(meas_theta))**2 - \
        (beam_maj*pl.sin(beam_theta))**2 - (beam_min*pl.cos(beam_theta))**2
  
    gamma = \
        2.*((meas_min**2-meas_maj**2)*pl.sin(meas_theta)*pl.cos(meas_theta) - 
            (beam_min**2-beam_maj**2)*pl.sin(beam_theta)*pl.cos(beam_theta))

    s = alpha + beta
    t = pl.sqrt((alpha-beta)**2 + gamma**2)

    #FIND THE SMALLEST RESOLUTION
    limit = pl.array([meas_min, meas_maj, beam_maj, beam_min]).min()
    limit = 0.1*limit*limit
  

    src_maj = pl.nan
    src_min = pl.nan
    src_pa = pl.nan
    
    if (alpha < 0 or beta < 0):
        # complete failure - alpha, beta are squares:
        worked = False

        #    ... CLOSE TO A POINT SOURCE
        if (0.5*(s-t) < limit and 
            alpha > -1*limit and 
            beta > -1*limit):
            point = True
            src_maj = 0.
            src_min = 0.
            src_pa = pl.nan
        else:
            point = False     

    else:
        src_maj = pl.sqrt(0.5*(s+t))
#        if (pl.absolute(gamma)+pl.absolute(alpha-beta) == 0):
#            src_pa = 0
#        else:
        src_pa  = 0.5*pl.arctan2(-1*gamma,alpha-beta)        

        if (s < t):
            # soft failure - de_maj may still be ok, but still call it failed
#            print alpha,beta,s,t
            worked=False
            point=False
            src_min=0.
        else:
            #... SUCCESS
            src_min = pl.sqrt(0.5*(s-t))
            worked = True
            point = False

    return src_maj,src_min,src_pa,point,worked



