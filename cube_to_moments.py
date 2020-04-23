#######################################################################
def cube_to_moments( datacube, assigncube, montecarlo=0, bm_pix=[0,0,0], verbose=False, fluxmap=None):
    
    mckeys=['mom0','ex_mom0','mom1v','mom1y','mom1x',
            'mom2v','ex_mom2v','mom2y','ex_mom2y','mom2x','ex_mom2x',
            'mom2maj','ex_mom2maj','mom2min','ex_mom2min',
            'max','halfmax_ell_maj','halfmax_ell_min','posang',
            'de_mom2maj','de_mom2min','de_posang',
            'de_halfmax_ell_maj','de_halfmax_ell_min',
            'max_integrated','fnu_maxchan','fnu_maxintchan','flux']
    
    # checks on inputs:
    if len(datacube.shape)!=3:
        print("data doesn't have correct shape - its ",datacube.shape)
    #    return
    if len(assigncube.shape)!=3:
        print("assign doesn't have correct shape - its ",assigncube.shape)
    #    return

    import matplotlib.pyplot as pl
    import numpy as np
    # create cloud list:
    zassigned=np.where(assigncube>0)
    if len(zassigned[0])==0:
        print("assignment cube appears to be blank")
    #    return
    cloudlist=np.unique(assigncube[zassigned])
    
    # vectorize the cubes for speed
    assign_vec = assigncube[zassigned]
    data_vec   = datacube[zassigned]
    # idl version used ind_to_xyz but python already has that as v=zassigned[0] etc
    # (remembering that cube is in v,y,x order in python)
    
    # to be used for montecarlo cube convolution
    beamsize_pix=np.sqrt(bm_pix[0]*bm_pix[1]) # convolution can only do circular for now
                
    def mad(data, axis=None):
        return np.nanmedian(np.absolute(data - np.nanmedian(data, axis)), axis)
    
    # measure the moments
    ncl=len(cloudlist)
    from measure_moments import measure_moments
    moments={}
    print("measuring moments")
    for i in range(ncl):
    #    if verbose: print "cloud  %5i / %5i \r"%(cloudlist[i],cloudlist.max()),
    
        zi=np.where( (assign_vec==cloudlist[i]) * (np.isnan(data_vec)==False) )[0]
        if len(zi)<=0:
            print("cloud ",cloudlist[i]," is bad")
        t=data_vec[zi].copy()
        v=zassigned[0][zi].copy()
        y=zassigned[1][zi].copy()
        x=zassigned[2][zi].copy()

        m=measure_moments(t,x,y,v,bm_pix=bm_pix)

        #import pdb
        #pl.clf()
        #pl.imshow(datacube[int(m['mom1v']),:,:],origin="bottom",interpolation="nearest")
        #pl.xlim(x.min(),x.max())
        #pl.ylim(y.min(),y.max())
        #pl.plot([m['mom1x']-m['mom2x'],m['mom1x']+m['mom2x']],[m['mom1y']-m['mom2y'],m['mom1y']+m['mom2y']])
        #pl.plot([m['mom1x']+m['halfmax_ell_maj'],m['mom1x']-m['halfmax_ell_maj']],[m['mom1y']-m['halfmax_ell_min'],m['mom1y']+m['halfmax_ell_min']])
        #pdb.set_trace()


        for k0,v0 in m.items():
            try:
                moments[k0].append(v0)
            except KeyError:
                moments[k0]=[v0]
        try:
            moments['id'].append(cloudlist[i])
        except KeyError:
            moments['id']=[cloudlist[i]]
            
        # gaussfit to spec:
        # reassemble spec at position of mom1x, mom1y
    #    mom1pos=[int(round(m['mom1x'])),int(round(m['mom1y']))]
    #    z=np.where( (x==mom1pos[0]) * (y==mom1pos[1]) )[0]
    #    spec=v[z]
    #    pdb.set_trace()
    
    # error on moments:
    mcmoments={}
    if montecarlo>0:
        import sys
        sys.stdout.write("noisy cube")
        from add_noise_to_cube import add_noise_to_cube
        mc=montecarlo
        for k in mckeys:
            mcmoments[k]=np.zeros([ncl,mc])
        for j in range(mc):
            sys.stdout.write("\rnoisy cube %03i/%03i"%(j,mc))
            sys.stdout.flush()
            noisy_cube=add_noise_to_cube(datacube,beamsize_pix,fluxmap=None)
            mcdata_vec = noisy_cube[zassigned]
            for i in range(ncl):
                zi=np.where( (assign_vec==cloudlist[i]) * (np.isnan(mcdata_vec)==False) )[0]
                v=zassigned[0][zi].copy()
                y=zassigned[1][zi].copy()
                x=zassigned[2][zi].copy()
                mj=measure_moments(mcdata_vec[zi].copy(),x,y,v,bm_pix=bm_pix)
                for k in mckeys:
                    mcmoments[k][i,j]=mj[k]
        print
        for i in range(ncl):
            for k in mckeys:
                try:
                    moments['d'+k][i]=mad(mcmoments[k][i])/.6745 # mad =.6745 sigma
                except KeyError:
                    moments['d'+k]=np.zeros(ncl)
                    moments['d'+k][i]=mad(mcmoments[k][i])/.6745 # mad =.6745 sigma
    
                # test:
                if mc>100:
                    try:
                        moments['d10'+k][i]=mad(mcmoments[k][i][0:10])/.6745 
                    except KeyError:
                        moments['d10'+k]=np.zeros(ncl)
                        moments['d10'+k][i]=mad(mcmoments[k][i][0:10])/.6745 
                    try:
                        moments['d30'+k][i]=mad(mcmoments[k][i][0:30])/.6745 
                    except KeyError:
                        moments['d30'+k]=np.zeros(ncl)
                        moments['d30'+k][i]=mad(mcmoments[k][i][0:30])/.6745 
                    try:
                        moments['d100'+k][i]=mad(mcmoments[k][i][0:100])/.6745 
                    except KeyError:
                        moments['d100'+k]=np.zeros(ncl)
                        moments['d100'+k][i]=mad(mcmoments[k][i][0:100])/.6745 
    
    for k,v in moments.items():
        moments[k]=np.array(v)
    
    return moments, mcmoments












