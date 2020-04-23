def pyprops(datafile,fluxfile,assignfile,root,assignfile2=None,montecarlo=0,doplots=False):

    import sys
    gitpaths=['/Users/remy/local/github/pyprops/']
    for gitpath in gitpaths:
        if not gitpath in sys.path:
            sys.path.insert(0,gitpath)
    
    print("importing modules and data")
    import astropy.io.fits as fits
    from astropy.table import Table
    
    from measure_moments import measure_moments
    from extrap import extrap
    from deconvolve_gauss import deconvolve_gauss
    from ellfit import ellfit
    from add_noise_to_cube import add_noise_to_cube
    
    import pylab as pl
    import numpy as np
    pl.ion()
    
    datacube=fits.getdata(datafile)
    s=datacube.shape
    if len(s)==4:
        if s[0]==1:
            #  cube has the trailing stokes axis which turns into a leading degenerate axis in python
            datacube=datacube[0]
    
    fluxmap=fits.getdata(fluxfile)
    
    
    hdr=fits.getheader(datafile)
    bmaj=hdr['bmaj'] # degrees
    bmin=hdr['bmin'] # degrees
    bpa=hdr['bpa']*np.pi/180 # ->rad
    
    from astropy import wcs
    w = wcs.WCS(hdr)
#    pix=pl.sqrt(-pl.det(w.celestial.pixel_scale_matrix)) # degrees
    pix=np.absolute(w.wcs.get_cdelt()[0]) # degrees
    
    bmaj_pix=bmaj/pix
    bmin_pix=bmin/pix
    bm_pix=[bmaj_pix,bmin_pix,bpa]
    
    import time
    x=time.localtime()
    datestr=str(x.tm_year)+("%02i"%x.tm_mon)+("%02i"%x.tm_mday)
    import pickle
    
    from cube_to_moments import cube_to_moments
    
    
    assigncube=fits.getdata(assignfile)
    moments,mcmoments=cube_to_moments(datacube,assigncube,montecarlo=montecarlo,bm_pix=bm_pix,fluxmap=fluxmap)
    
    moments['posang']=moments['posang']%np.pi
    moments['de_posang']=moments['de_posang']%np.pi
    
    pickle.dump([moments,bm_pix,mcmoments],open(root+".pyprops.pkl","wb"))

    
    if assignfile2:
        assigncube2=fits.getdata(assignfile2)
        moments2,mcmoments2=cube_to_moments(datacube,assigncube2,montecarlo=montecarlo,bm_pix=bm_pix,fluxmap=fluxmap)
        
        moments2['posang']=moments2['posang']%np.pi
        moments2['de_posang']=moments2['de_posang']%np.pi
        
        pickle.dump([moments2,bm_pix,mcmoments2],open(root2+"pyprops.pkl","wb"))


    #######################################################################
    # diagnostic plots

    if doplots:
        pl.clf()
        ellrad=np.sqrt(moments['halfmax_ell_maj']*moments['halfmax_ell_min'])
        pl.plot(ellrad,moments['mom2v'],'.')
        
        if assignfile2:
            ellrad2=np.sqrt(moments2['halfmax_ell_maj']*moments2['halfmax_ell_min'])
            pl.plot(ellrad2,moments2['mom2v'],'.')
        
        pl.xlabel("size")
        pl.ylabel("linewidth")
        
        pl.xscale("log")
        pl.yscale("log")
        
        pl.savefig(root+".pyprops.sizeline.png")


    if doplots:
        pl.clf()
        u=np.argsort(moments['flux'])[::-1]
        pl.plot((moments['flux'][u]),pl.arange(len(u)))
        if assignfile2:
            u2=np.argsort(moments2['flux'])[::-1]
            pl.plot((moments2['flux'][u2]),np.arange(len(u2)))
        pl.yscale("log")
        pl.xscale("log")
        pl.xlabel("flux")
        pl.xlim(pl.xlim()[::-1])
        pl.ylabel("cumulative number")
        
        pl.savefig(root+".pyprops.cumfluxdist.png")
        

    if doplots and montecarlo>0:
        pl.clf()
        fnu=0.5*(moments['fnu_maxintchan']+moments['fnu_maxchan'])
        delfnu=np.absolute(moments['fnu_maxintchan']-moments['fnu_maxchan'])
        dfnu=0.5*np.sqrt( moments['dfnu_maxchan']**2 + moments['dfnu_maxintchan']**2 )
        z=np.where(delfnu>dfnu)[0]
        if len(z)>0:
            dfnu[z]=delfnu[z]
        
        pl.errorbar(fnu,moments['flux'],xerr=dfnu,yerr=moments['dflux'],fmt='.')
        pl.xlabel("I [Bunit * vpix]")
        pl.ylabel("F [Bunit * pix^2 * vpix")
        pl.xscale("log")
        pl.yscale("log")
        pl.xlim(.1,100)
        pl.ylim(.5,500)
        pl.plot(pl.xlim(),pl.array(pl.xlim())*3)
        pl.savefig(root+".pyprops.I_F.png")




    # size comparisons
    if doplots and montecarlo>0:
        pl.clf()
        bmarea=bmaj_pix*bmin_pix*np.pi/4 # not a beam "volume" 
        area1=moments['de_mom2maj']*moments['de_mom2min']*2.354**2*np.pi/4/bmarea
        darea1=area1*np.sqrt( (moments['dde_mom2maj']/moments['de_mom2maj'])**2+
                              (moments['dde_mom2min']/moments['de_mom2min'])**2 )
        area2=moments['de_halfmax_ell_maj']*moments['de_halfmax_ell_min']*np.pi/4/bmarea
        darea2=area2*np.sqrt( (moments['dde_halfmax_ell_maj']/moments['de_halfmax_ell_maj'])**2+
                              (moments['dde_halfmax_ell_min']/moments['de_halfmax_ell_min'])**2 )
        
        pl.errorbar(moments['flux'],area1,xerr=moments['dflux'],yerr=darea1,fmt='.',label="mom2")
        pl.errorbar(moments['flux'],area2,xerr=moments['dflux'],yerr=darea2,fmt='.',label="ellfit")
        pl.ylabel("area [beams]")
        pl.xlabel("F [Bunit * pix^2 * vpix]")
        pl.legend(loc="best",prop={"size":10},numpoints=1)
        pl.xscale("log")
        pl.savefig(root+".pyprops.F_area.png")
        
        
        
        bmarea=bmaj_pix*bmin_pix*np.pi/4 # not a beam "volume" 
        area1=moments['mom2maj']*moments['mom2min']*2.354**2*np.pi/4/bmarea
        darea1=area1*np.sqrt( (moments['dmom2maj']/moments['mom2maj'])**2+
                              (moments['dmom2min']/moments['mom2min'])**2 )
        de_area1=moments['de_mom2maj']*moments['de_mom2min']*2.354**2*np.pi/4/bmarea
        dde_area1=de_area1*np.sqrt( (moments['dde_mom2maj']/moments['de_mom2maj'])**2+
                              (moments['dde_mom2min']/moments['de_mom2min'])**2 )
        area2=moments['halfmax_ell_maj']*moments['halfmax_ell_min']*np.pi/4/bmarea
        darea2=area2*np.sqrt( (moments['dhalfmax_ell_maj']/moments['halfmax_ell_maj'])**2+
                              (moments['dhalfmax_ell_min']/moments['halfmax_ell_min'])**2 )
        de_area2=moments['de_halfmax_ell_maj']*moments['de_halfmax_ell_min']*np.pi/4/bmarea
        dde_area2=de_area2*np.sqrt( (moments['dde_halfmax_ell_maj']/moments['de_halfmax_ell_maj'])**2+
                              (moments['dde_halfmax_ell_min']/moments['de_halfmax_ell_min'])**2 )
        
        pl.clf()
        pl.subplot(211)
        pl.errorbar(area1,de_area1,xerr=darea1,yerr=dde_area1,fmt='.',label="success")
        z=np.where(np.isnan(de_area1))[0]
        pl.errorbar(area1[z],area1[z],xerr=darea1[z],fmt='.',label="failed")
        z=np.where(de_area1==0)[0]
        pl.errorbar(area1[z],area1[z],xerr=darea1[z],fmt='.',label="ptsrc")
        pl.ylabel("area [beams, deconv]")
        pl.xlabel("area [beams, meas]")
        pl.legend(loc="best",prop={"size":10},numpoints=1)
        
        pl.subplot(212)
        pl.errorbar(area2,de_area2,xerr=darea2,yerr=dde_area2,fmt='.',label="success")
        z=np.where(np.isnan(de_area2))[0]
        pl.errorbar(area2[z],area2[z],xerr=darea2[z],fmt='.',label="failed")
        z=np.where(de_area2==0)[0]
        pl.errorbar(area2[z],area2[z],xerr=darea2[z],fmt='.',label="ptsrc")
        pl.ylabel("area [beams, deconv]")
        pl.xlabel("area [beams, meas]")
        pl.legend(loc="best",prop={"size":10},numpoints=1)
        
        #pl.xscale("log")
        pl.savefig(root+".pyprops.area_dearea.png")
        
        
        raterr=area2/area1
        draterr=raterr*np.sqrt( (darea2/area2)**2 + (darea1/area1)**2 )
        pl.clf()
        pl.errorbar(area2,raterr,xerr=darea2,yerr=draterr,fmt='.')
        pl.xlabel("ell area")
        pl.ylabel("ell area / mom2 area")
    
        pl.savefig(root+".pyprops.ellarea_momarea.png")
    
    
    
    # fluxes and integrated spectra:
    if doplots and montecarlo>0:
        ratiofnu=moments['fnu_maxintchan']/moments['fnu_maxchan']
        dratio=ratiofnu*np.sqrt( (moments['dfnu_maxchan']/moments['fnu_maxchan'])**2 +
                                 (moments['dfnu_maxintchan']/moments['fnu_maxintchan'])**2 )
        difffnu=moments['fnu_maxintchan']-moments['fnu_maxchan']
        ddiff=np.sqrt( moments['dfnu_maxchan']**2 +
                       moments['dfnu_maxintchan']**2 )
        
        fnu=0.5*(moments['fnu_maxintchan']+moments['fnu_maxchan'])
        
        pl.clf()
        z=np.where(difffnu>0)[0]
        pl.errorbar(moments['fnu_maxchan'][z],(difffnu/fnu)[z],xerr=moments['dfnu_maxchan'][z],yerr=(ddiff/fnu)[z],fmt='.')
        pl.xlabel("Fnu @max    [Bunit * pix^2]")
        pl.ylabel("(Fnu @maxint - Fnu @max)/Fnu")
        pl.xscale("log")
        pl.xlim(.3,40)
        pl.ylim(-.3,1.5)
        pl.savefig("pyprops."+datestr+".F_maxint_v_max.png")
    
    
    # convergence of errors with MC iterations:
    if doplots and montecarlo>100:
        pl.clf()
        k='max'
        k='de_mom2min'
        if montecarlo>100:
            pl.plot(moments[k],moments['d10'+k],'.',label='10')
            pl.plot(moments[k],moments['d30'+k],'.',label='30')
            pl.plot(moments[k],moments['d100'+k],'.',label='100')
        pl.plot(moments[k],moments['d'+k],'.',label=str(montecarlo))
        pl.legend(loc="best",prop={"size":10},numpoints=1)
        pl.xlabel(k)
        pl.ylabel("d"+k)
        
        if k=='max':
            pl.plot(pl.xlim(),[0.003,0.003],'k')
            pl.plot(pl.xlim(),[0.005,0.005],'k',linestyle="dashed")
            pl.xlim(0,0.35)
        
        if montecarlo>100:
            pl.savefig(root+".pyprops."+k+".mciters.png")
        else:
            pl.savefig(root+".pyprops."+k+".png")        


    if doplots:
        # previously, we did
        # if maj<bmaj*1.1: set de_maj,de_min to bmaj/2, bmin/2
        # elif min<bmin*1.1: de_maj=pl.sqrt(maj**2-bmaj**2), de_min=bmin/2
        
        pl.clf()
        pl.plot(moments['halfmax_ell_maj'],moments['halfmax_ell_min'],'.')
        pl.xlabel("measured major ell @halfmax")
        pl.ylabel("measured minor ell @halfmax")
        pl.plot([bm_pix[0],pl.xlim()[1]],[bm_pix[1],bm_pix[1]],'k',linestyle='dotted')
        pl.plot([bm_pix[0],bm_pix[0]],[bm_pix[1],pl.ylim()[1]],'k',linestyle='dotted')
        # if its something half the beamsize, convolved with the beam, 
        # it'll end up with size=pl.sqrt(1+.5**2)=1.12*bm
        pl.plot([bm_pix[0]*1.12,pl.xlim()[1]],[bm_pix[1]*1.12,bm_pix[1]*1.12],'k',linestyle='dotted')
        pl.plot([bm_pix[0]*1.12,bm_pix[0]*1.12],[bm_pix[1]*1.12,pl.ylim()[1]],'k',linestyle='dotted')
        
        z=np.where(np.isnan(moments['de_halfmax_ell_maj']))[0]
        pl.plot(moments['halfmax_ell_maj'][z],moments['halfmax_ell_min'][z],'s',label='dec=nan')
        z=np.where((moments['de_halfmax_ell_min']==0)*(moments['de_halfmax_ell_maj']>0))[0]
        pl.plot(moments['halfmax_ell_maj'][z],moments['halfmax_ell_min'][z],'cd',label='dec min=0')
        z=np.where((moments['de_halfmax_ell_min']==0)*(moments['de_halfmax_ell_maj']==0))[0]
        pl.plot(moments['halfmax_ell_maj'][z],moments['halfmax_ell_min'][z],'r*',label='dec both=0')
        pl.legend(loc="best",prop={"size":10},numpoints=1)
        
        pl.savefig(root+".pyprops.measured_ellipses.png")
        
        
        pl.clf()
        pl.plot(moments['posang'],moments['halfmax_ell_min'],'.')
        pl.xlabel("posang")
        pl.ylabel("measured minor ell @halfmax")
        pl.plot([bm_pix[2],bm_pix[2]],pl.ylim(),'k',linestyle='dashed')
        pl.plot(pl.xlim(),[bm_pix[1],bm_pix[1]],'k',linestyle='dotted')
        pl.plot(pl.xlim(),[bm_pix[1]*1.12,bm_pix[1]*1.12],'k',linestyle='dotted')
        
        z=np.where(np.isnan(moments['de_halfmax_ell_maj']))[0]
        pl.plot(moments['posang'][z],moments['halfmax_ell_min'][z],'s',label='dec=nan')
        z=np.where((moments['de_halfmax_ell_min']==0)*(moments['de_halfmax_ell_maj']>0))[0]
        pl.plot(moments['posang'][z],moments['halfmax_ell_min'][z],'cd',label='dec min=0')
        z=np.where((moments['de_halfmax_ell_min']==0)*(moments['de_halfmax_ell_maj']==0))[0]
        pl.plot(moments['posang'][z],moments['halfmax_ell_min'][z],'r*',label='dec both=0')
        pl.legend(loc="best",prop={"size":10},numpoints=1)
        
        pl.savefig(root+".pyprops.measured_ellipses_angles.png")
        
        
        
        p  =moments['posang'].copy()
        dp=moments['de_posang'].copy()
        p=(p-bm_pix[2]+np.pi/2)%np.pi +bm_pix[2]-np.pi/2
        dp=(dp-bm_pix[2]+np.pi/2)%np.pi +bm_pix[2]-np.pi/2
        
        pl.clf()
        pl.plot(p,dp,'.')
        pl.plot([bm_pix[2],bm_pix[2]],pl.ylim(),'k')
        pl.plot(pl.xlim(),[bm_pix[2],bm_pix[2]],'k')
        pl.xlabel("meas posang")
        pl.ylabel("deconv posang")
        
        pl.plot(pl.xlim(),[bm_pix[2]+np.pi/2,bm_pix[2]+np.pi/2],'k',linestyle='dotted')
    
        pl.savefig(root+".pyprops.dec_ellipses_angles.png")
