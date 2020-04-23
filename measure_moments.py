def measure_moments(tin,xin,yin,vin, bm_pix=None , do_extrap=True, verbose=False, maxindex=None):
    #import pdb
    # for testing:
    if not maxindex: maxindex=len(tin)
    import matplotlib.pyplot as pl
    import numpy as np
    
    # sort by pixel values from brightest to faintest
    u=np.argsort(tin)[::-1][0:maxindex-1]
    t=tin[u].copy()
    v=vin[u].copy()  # these are all cube indices in pixels
    y=yin[u].copy()
    x=xin[u].copy()

    if len(t)==0:
        t=np.array([-1])
        v=np.array([-1])
        y=np.array([-1])
        x=np.array([-1])
    
    moment={}
    moment['npix']=len(t)
    moment['max']=t.max()
    moment['min']=t.min()
    moment['posmax']=(v[0],y[0],x[0]) # respect pythonic order of v,y,x so this can be used directly to index in the cube 
    
    # TODO IMPLEMENT DO_CLIP (remove edge value of cloud)
    
    # use cumulative moments (need a function that returns the cumulative array)
    # or else just return the sum
    if do_extrap:
        tot=np.cumsum  
    else:
        tot=np.sum   # TODO may break some things below with indexing mom1x
    
    if t.max()>0:
    
        # moments
        mom0  = tot(t) 
        mom1v = tot(v*t) / mom0
        mom1y = tot(y*t) / mom0
        mom1x = tot(x*t) / mom0
        
        # second moments
        term1x = tot(t*x**2)
        term2x = tot(t*x)**2 / mom0
        mom2x = np.sqrt((term1x - term2x)/mom0)
        zeroind = np.where(np.absolute(term1x - term2x) < 1e-10)
        if len(zeroind)>0: mom2x[zeroind] = 0.0
        
        term1y = tot(t*y**2)
        term2y = tot(t*y)**2 / mom0
        mom2y = np.sqrt((term1y - term2y)/mom0)
        zeroind = np.where(np.absolute(term1y - term2y) < 1e-10)
        if len(zeroind)>0: mom2y[zeroind] = 0.0
        
        term1v = tot(t*v**2)
        term2v = tot(t*v)**2 / mom0
        mom2v = np.sqrt((term1v - term2v)/mom0)
        zeroind = np.where(np.absolute(term1v - term2v) < 1e-10)
        if len(zeroind)>0: mom2v[zeroind] = 0.0
          
        
        ###########################################################
        # area, vel width, and ellipse fit
        
        # unique set of 2-d indices for the 3-d clump
        twod_id = x + y*(x.max()+1)  
        twod_ind = np.unique(twod_id,return_index=True)[1]
        twod_x = x[twod_ind]
        twod_y = y[twod_ind]
        
        area=len(twod_ind)

        from ellfit import ellfit
        # major/minor axis fit to full projected extent of the clump:
        ell_maj,ell_min,ell_pa=ellfit(twod_x,twod_y)  
        # this really is a major/minor axis (factor of 4 in my new ellfit.py)
        
        # now fit the entire 3d clump, but weighted, to get a position angle:
        well_maj,well_min,posang=ellfit(x,y,wt=t)
        # the major and minor axes here are related to the second moments, 
        # but only if the assignment region goes down to near-zero intensity.
        # however, the position angle is robust.  
        
        
        #-----------------------
        # now find all pixels above 0.5*max, and again find the 2D projection:
        half_ind = np.where(t>0.5*t.max())[0]
        half_twod_ind = np.unique(twod_id[half_ind],return_index=True)[1]
        half_twod_x = x[half_twod_ind]
        half_twod_y = y[half_twod_ind]
        
        half_area=len(half_twod_ind)

        half_ell_maj,half_ell_min,half_ell_pa=ellfit(half_twod_x,half_twod_y)  
                
        
        # next, use the 3d-calculated posang to fit the moments along the major/minor
        xrot =  x*np.cos(posang) + y*np.sin(posang)
        yrot = -x*np.sin(posang) + y*np.cos(posang)    
        
        mom1yrot = tot(yrot*t) / mom0
        mom1xrot = tot(xrot*t) / mom0
        
        # second moments
        term1xrot = tot(t*xrot**2)
        term2xrot = tot(t*xrot)**2 / mom0
        mom2xrot = np.sqrt((term1xrot - term2xrot)/mom0)
        zeroind = np.where(np.absolute(term1xrot - term2xrot) < 1e-10)
        if len(zeroind)>0: mom2xrot[zeroind] = 0.0
        
        term1yrot = tot(t*yrot**2)
        term2yrot = tot(t*yrot)**2 / mom0
        mom2yrot = np.sqrt((term1yrot - term2yrot)/mom0)
        zeroind = np.where(np.absolute(term1yrot - term2yrot) < 1e-10)
        if len(zeroind)>0: mom2yrot[zeroind] = 0.0
    



    #-----------------
    # Inu=max
    # I  = int(Inu dnu)
    # max is at posmax (=x[0],y[0])
    z=np.where( (x==x[0])*(y==y[0]) )[0]
    # vpix=vin[z]
    # vspec=tin[z]
    moment['max_integrated']=t[z].sum()  # (input units)*chanwid

    # Fnu = int(Inu*dA)  at the peak channel, but what chan?
    maxchan=v[0]
    z=np.where( v==v[0] )[0]
    moment['fnu_maxchan']=t[z].sum() # (input units)*pixels

    # integrated spectrum
    vpix=np.unique(v)
    vpix.sort()
    nvpix=len(vpix)
    intspec=np.zeros(nvpix)
    for i in range(nvpix):
        z=np.where( v==vpix[i] )[0]
        intspec[i]=t[z].sum()
    maxintchan=vpix[np.where( intspec==intspec.max() )[0]]
    z=np.where( v==maxintchan )[0]
    moment['fnu_maxintchan']=t[z].sum()

    # F   = int(Inu*dA*dnu) = total
    moment['flux']=tin.sum()

    #######################################
    moment['vrange'] = (v.min(),v.max())
    
    if t.max()>0: 
        moment['assign_area'] = area
        moment['mom0']  =  mom0[-1]
        moment['mom1x'] = mom1x[-1]
        moment['mom1y'] = mom1y[-1]
        moment['mom1v'] = mom1v[-1]
        moment['mom2x'] = mom2x[-1]
        moment['mom2y'] = mom2y[-1]
        moment['mom2v'] = mom2v[-1]
    else:
        moment['assign_area'] = 0
        moment['mom0']  = 0.
        moment['mom1x'] = 0.
        moment['mom1y'] = 0.
        moment['mom1v'] = 0.
        moment['mom2x'] = 0.
        moment['mom2y'] = 0.
        moment['mom2v'] = 0.
        
    if t.max()>0:
        # major/minor axes of assign projected 2d area
        moment['assign_ell_maj'] = ell_maj 
        moment['assign_ell_min'] = ell_min
    
        # pos angle of weighted all points in 3d - used for major/minor moments
        moment['posang'] = posang # radians
        
        moment['halfmax_ell_maj'] = half_ell_maj
        moment['halfmax_ell_min'] = half_ell_min
        
        moment['halfmax_area'] = half_area
        moment['halfmax_vrange'] = (v[half_ind].min(),v[half_ind].max())
        
        if mom2xrot[-1]>mom2yrot[-1]:
            moment['mom2maj'] = mom2xrot[-1]
            moment['mom2min'] = mom2yrot[-1]
        else:                           
            moment['mom2maj'] = mom2yrot[-1]
            moment['mom2min'] = mom2xrot[-1]

    else:
        moment['assign_ell_maj'] = 0.
        moment['assign_ell_min'] = 0.
        moment['posang'] = 0.        
        posang = 0. # for deconvolution below
        half_ell_maj=0.
        half_ell_min=0.
        moment['halfmax_ell_maj'] = 0.
        moment['halfmax_ell_min'] = 0.
        moment['halfmax_area'] = 0.
        moment['halfmax_vrange'] = (-1,-1)
        moment['mom2maj'] = 0.
        moment['mom2min'] = 0.
    
    
    ###############################################
    if do_extrap and t.max()>0:
        from extrap import extrap
        ex_mom2x = extrap(t,mom2x)
        ex_mom2y = extrap(t,mom2y)
        ex_mom2v = extrap(t,mom2v)
        ex_mom2xrot = extrap(t,mom2xrot)
        ex_mom2yrot = extrap(t,mom2yrot)
        ex_mom0 = extrap(t,mom0)

        # TODO extrapolated posang ?
    
        if ex_mom0<np.nanmax(mom0):
#            print "ERROR in EXTRAPOLATION!!! Extrapolated flux is LESS THAN measured flux. That is BAD!"
#            print 'Defaulting to unweighted extrapolation...'
            ex_mom0 = extrap(t,mom0,weight=False)
        
        moment['ex_mom0']     = ex_mom0 
        moment['ex_mom2x']    = ex_mom2x 
        moment['ex_mom2y']    = ex_mom2y 
        moment['ex_mom2v']    = ex_mom2v 
        if mom2xrot[-1]>mom2yrot[-1]:
            moment['ex_mom2maj'] = ex_mom2xrot
            moment['ex_mom2min'] = ex_mom2yrot
        else:
            moment['ex_mom2maj'] = ex_mom2yrot
            moment['ex_mom2min'] = ex_mom2xrot

    else:
        moment['ex_mom0']     = 0.
        moment['ex_mom2x']    = 0.
        moment['ex_mom2y']    = 0.
        moment['ex_mom2v']    = 0.
        moment['ex_mom2maj']  = 0.
        moment['ex_mom2min']  = 0.





    #######################################
    # deconolution
    if bm_pix:
        if do_extrap:
            ma=moment['ex_mom2maj']
            mi=moment['ex_mom2min']
        else:
            ma=moment['mom2maj']
            mi=moment['mom2min']

        from deconvolve_gauss import deconvolve_gauss
        de_mom2maj,de_mom2min,de_posang,pt,worked = deconvolve_gauss(ma*2.354,mi*2.354,posang,bm_pix[0],bm_pix[1],bm_pix[2])
        moment['de_mom2maj']=de_mom2maj/2.354
        moment['de_mom2min']=de_mom2min/2.354
        moment['de_posang']=de_posang

        de_hellmaj,de_hellmin,de_hellposang,pt,worked = deconvolve_gauss(half_ell_maj,half_ell_min,posang,bm_pix[0],bm_pix[1],bm_pix[2])
        moment['de_halfmax_ell_maj']=de_hellmaj
        moment['de_halfmax_ell_min']=de_hellmin

    elif verbose:
        print("can't do deconvolved quantities - beam not known")

        
    
    if verbose:
        print("0  :", moment['mom0'],moment['ex_mom0'])
        print("2x :", moment['mom2x'],moment['ex_mom2x'])
        print("2y :", moment['mom2y'],moment['ex_mom2y'])
        print("2v :", moment['mom2v'],moment['ex_mom2v'])
        print("2ma:", moment['mom2maj'],moment['ex_mom2maj'])
        print("2mi:", moment['mom2min'],moment['ex_mom2min'])
    
    return moment
