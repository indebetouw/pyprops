dir='/Users/remy/lustre/naasc/users/rindebet/30dor/cycle0/'
datafile=dir+'12CO.combined.20130113.nopbcor.Jybm.kms.fits'
fluxfile=dir+'cal.contsub.02.lineself.13CO.cb.1st.p075.flux.ch30.fits'
assignfile=dir+'20130410.assign_out.fits'
root="dor.c0.12co.2013assign"

assignfile2=dir+'12CO.combined.20130113.nopbcor.Jybm.kms.clumps.fits'
root2="dor.c0.12co.2018assign"


montecarlo=0 # 300
doplots=False

# force recalc:
# del moments

import sys
gitpaths=['/Users/remy/lustre/naasc/users/rindebet/github/pyprops/']
for gitpath in gitpaths:
    if not gitpath in sys.path:
        sys.path.insert(0,gitpath)

print "importing modules and data"
import astropy.io.fits as fits
from astropy.table import Table

from measure_moments import measure_moments
from extrap import extrap
from deconvolve_gauss import deconvolve_gauss
from ellfit import ellfit
from add_noise_to_cube import add_noise_to_cube

import pylab as pl
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
bpa=hdr['bpa']*pl.pi/180 # ->rad

from astropy import wcs
w = wcs.WCS(hdr)
pix=pl.sqrt(-pl.det(w.celestial.pixel_scale_matrix)) # degrees

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

moments['posang']=moments['posang']%pl.pi
moments['de_posang']=moments['de_posang']%pl.pi

pickle.dump([moments,bm_pix,mcmoments],open("pyprops."+datestr+"."+root+".pkl","wb"))


assigncube2=fits.getdata(assignfile2)
moments2,mcmoments2=cube_to_moments(datacube,assigncube2,montecarlo=montecarlo,bm_pix=bm_pix,fluxmap=fluxmap)

moments2['posang']=moments2['posang']%pl.pi
moments2['de_posang']=moments2['de_posang']%pl.pi

pickle.dump([moments2,bm_pix,mcmoments2],open("pyprops."+datestr+"."+root2+".pkl","wb"))



pl.clf()
ellrad=pl.sqrt(moments['halfmax_ell_maj']*moments['halfmax_ell_min'])
pl.plot(ellrad,moments['mom2v'],'.')

ellrad2=pl.sqrt(moments2['halfmax_ell_maj']*moments2['halfmax_ell_min'])
pl.plot(ellrad2,moments2['mom2v'],'.')

pl.xlabel("size")
pl.ylabel("linewidth")

pl.xscale("log")
pl.yscale("log")

pl.savefig("pyprops."+datestr+".sizeline.2assign.png")



pl.clf()
#h1=pl.histogram(pl.log10(moments['flux'] ),bins=0.8+pl.frange(10)*0.45)
#h2=pl.histogram(pl.log10(moments2['flux']),bins=0.8+pl.frange(10)*0.45)
#pl.plot(h1[1][1:][::-1],h1[0][::-1].cumsum())
#pl.plot(h2[1][1:][::-1],h2[0][::-1].cumsum())

pl.clf()
u=pl.argsort(moments['flux'])[::-1]
pl.plot((moments['flux'][u]),pl.arange(len(u)))
u2=pl.argsort(moments2['flux'])[::-1]
pl.plot((moments2['flux'][u2]),pl.arange(len(u2)))
pl.yscale("log")
pl.xscale("log")
pl.xlabel("flux")
pl.xlim(pl.xlim()[::-1])
pl.ylabel("cumulative number")


pl.savefig("pyprops."+datestr+".cumfluxdist.2assign.png")







#######################################################################
# diagnostic plots
if doplots and montecarlo>0:
    pl.clf()
    fnu=0.5*(moments['fnu_maxintchan']+moments['fnu_maxchan'])
    delfnu=pl.absolute(moments['fnu_maxintchan']-moments['fnu_maxchan'])
    dfnu=0.5*pl.sqrt( moments['dfnu_maxchan']**2 + moments['dfnu_maxintchan']**2 )
    z=pl.where(delfnu>dfnu)[0]
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
    pl.savefig("pyprops."+datestr+"."+root+".I_F.png")




# size comparisons
if doplots and montecarlo>0:
    pl.clf()
    bmarea=bmaj_pix*bmin_pix*pl.pi/4 # not a beam "volume" 
    area1=moments['de_mom2maj']*moments['de_mom2min']*2.354**2*pl.pi/4/bmarea
    darea1=area1*pl.sqrt( (moments['dde_mom2maj']/moments['de_mom2maj'])**2+
                          (moments['dde_mom2min']/moments['de_mom2min'])**2 )
    area2=moments['de_halfmax_ell_maj']*moments['de_halfmax_ell_min']*pl.pi/4/bmarea
    darea2=area2*pl.sqrt( (moments['dde_halfmax_ell_maj']/moments['de_halfmax_ell_maj'])**2+
                          (moments['dde_halfmax_ell_min']/moments['de_halfmax_ell_min'])**2 )
    
    pl.errorbar(moments['flux'],area1,xerr=moments['dflux'],yerr=darea1,fmt='.',label="mom2")
    pl.errorbar(moments['flux'],area2,xerr=moments['dflux'],yerr=darea2,fmt='.',label="ellfit")
    pl.ylabel("area [beams]")
    pl.xlabel("F [Bunit * pix^2 * vpix]")
    pl.legend(loc="best",prop={"size":10},numpoints=1)
    pl.xscale("log")
    pl.savefig("pyprops."+datestr+"."+root+".F_area.png")
    
    
    
    bmarea=bmaj_pix*bmin_pix*pl.pi/4 # not a beam "volume" 
    area1=moments['mom2maj']*moments['mom2min']*2.354**2*pl.pi/4/bmarea
    darea1=area1*pl.sqrt( (moments['dmom2maj']/moments['mom2maj'])**2+
                          (moments['dmom2min']/moments['mom2min'])**2 )
    de_area1=moments['de_mom2maj']*moments['de_mom2min']*2.354**2*pl.pi/4/bmarea
    dde_area1=de_area1*pl.sqrt( (moments['dde_mom2maj']/moments['de_mom2maj'])**2+
                          (moments['dde_mom2min']/moments['de_mom2min'])**2 )
    area2=moments['halfmax_ell_maj']*moments['halfmax_ell_min']*pl.pi/4/bmarea
    darea2=area2*pl.sqrt( (moments['dhalfmax_ell_maj']/moments['halfmax_ell_maj'])**2+
                          (moments['dhalfmax_ell_min']/moments['halfmax_ell_min'])**2 )
    de_area2=moments['de_halfmax_ell_maj']*moments['de_halfmax_ell_min']*pl.pi/4/bmarea
    dde_area2=de_area2*pl.sqrt( (moments['dde_halfmax_ell_maj']/moments['de_halfmax_ell_maj'])**2+
                          (moments['dde_halfmax_ell_min']/moments['de_halfmax_ell_min'])**2 )
    
    pl.clf()
    pl.subplot(211)
    pl.errorbar(area1,de_area1,xerr=darea1,yerr=dde_area1,fmt='.',label="success")
    z=pl.where(pl.isnan(de_area1))[0]
    pl.errorbar(area1[z],area1[z],xerr=darea1[z],fmt='.',label="failed")
    z=pl.where(de_area1==0)[0]
    pl.errorbar(area1[z],area1[z],xerr=darea1[z],fmt='.',label="ptsrc")
    pl.ylabel("area [beams, deconv]")
    pl.xlabel("area [beams, meas]")
    pl.legend(loc="best",prop={"size":10},numpoints=1)
    
    pl.subplot(212)
    pl.errorbar(area2,de_area2,xerr=darea2,yerr=dde_area2,fmt='.',label="success")
    z=pl.where(pl.isnan(de_area2))[0]
    pl.errorbar(area2[z],area2[z],xerr=darea2[z],fmt='.',label="failed")
    z=pl.where(de_area2==0)[0]
    pl.errorbar(area2[z],area2[z],xerr=darea2[z],fmt='.',label="ptsrc")
    pl.ylabel("area [beams, deconv]")
    pl.xlabel("area [beams, meas]")
    pl.legend(loc="best",prop={"size":10},numpoints=1)
    
    #pl.xscale("log")
    pl.savefig("pyprops."+datestr+"."+root+".area_dearea.png")
    
    
    raterr=area2/area1
    draterr=raterr*pl.sqrt( (darea2/area2)**2 + (darea1/area1)**2 )
    pl.clf()
    pl.errorbar(area2,raterr,xerr=darea2,yerr=draterr,fmt='.')
    pl.xlabel("ell area")
    pl.ylabel("ell area / mom2 area")

    pl.savefig("pyprops."+datestr+"."+root+".ellarea_momarea.png")



# fluxes and integrated spectra:
if doplots and montecarlo>0:
    ratiofnu=moments['fnu_maxintchan']/moments['fnu_maxchan']
    dratio=ratiofnu*pl.sqrt( (moments['dfnu_maxchan']/moments['fnu_maxchan'])**2 +
                             (moments['dfnu_maxintchan']/moments['fnu_maxintchan'])**2 )
    difffnu=moments['fnu_maxintchan']-moments['fnu_maxchan']
    ddiff=pl.sqrt( moments['dfnu_maxchan']**2 +
                   moments['dfnu_maxintchan']**2 )
    
    fnu=0.5*(moments['fnu_maxintchan']+moments['fnu_maxchan'])
    
    pl.clf()
    z=pl.where(difffnu>0)[0]
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
        pl.savefig("pyprops."+datestr+"."+root+"."+k+".mciters.png")
    else:
        pl.savefig("pyprops."+datestr+"."+root+"."+k+".png")


if False:
    # compare with previous
    
    oldposang=[2.4314,2.7109,2.1405,2.92264,3.09558,2.38106,2.52803,2.55262,2.02719,
               0.1892,3.0071,2.5895,2.50012,1.12271,1.07479,2.08484,1.91065,2.27928,
               2.5875,2.5108,2.6399,2.32694,2.97419,2.53362,1.85333,1.91494,3.04721,
               2.6266,1.6787,2.7401,2.20165,2.82723,0.84754,2.61641,0.73846,1.88631,
               2.2509,1.8303,2.2596,2.03882,1.96760,1.72844,2.02982,1.97079,2.06588,
               3.0506,1.6413,1.1097,1.74695,2.46724,2.72892,2.34758,2.49345,2.06025,
               2.2601,2.7932,2.5720,3.11948,2.49226,2.20629,2.96764,2.45699,2.27071,
               2.7325,2.4506,0.4373,2.08349,1.25807,2.51510,2.64043,1.79829,2.60664,
               2.5151,1.8262,2.8365,2.29365,2.83030,2.42320,2.45566,1.73173,2.14717,
               2.0150,1.7740,2.4085,2.44089,1.76299,2.70881,2.02988,3.07241,2.28563,
               1.5299,2.3213,0.8095,2.56525,2.16541,1.64124,1.84096,2.47796,0.29105,
               0.5811,2.0552,0.1292,2.38877,0.50973,1.63598,2.90378,2.72030,1.81256,
               1.9005,1.9151,1.0555,2.40282,1.71604,0.09452,0.05322,3.13331,0.11797,1.71203]
    
    t=Table.read("hcop.20160901.props3_2016.txt",format="ascii")
    # there are some missing ids from cloudlist, and more missing from t    
#    z0=pl.where(t['i']==31)[0]
#    oldposang=pl.append(oldposang[0:z0],oldposang[z0+1:])
#    z=pl.where(t['i']!=31)[0]
#    t=t[z]
    tindex=pl.zeros(len(t),dtype=int)-1
    for i in range(len(t)):
        z=pl.where(moments['id']==t['i'][i])[0]
        if len(z)>0:
            tindex[i]=z
    tfound=pl.where(tindex>=0)[0]
    tindex=tindex[tfound]
    
#    pl.clf()
#    pl.plot(pl.pi-oldposang,moments['posang'][tindex],'.')
    # bpa=0.88, most of oldposang are 2-2.5 i.e. pl.pi-pa ???
    # related to order of axes in ellfit.py arctan calculation
    
    
    
    pl.clf()
    #pl.plot(t['max'],moments['max'][tindex],'.') 
    k=['m_dv0','mom2v']
    pix=1.
    
    k=['m_maj0','mom2maj']
    k=['m_min0','mom2min']
    k=['m_maj2','de_mom2maj']
    k=['m_min2','de_mom2min']  # this one looks bad
    pix=206265./50000/.3 # from pc to arcsec and then to pix
    
    k=['e_maj0','halfmax_ell_maj']
    k=['e_min0','halfmax_ell_min']
    k=['e_maj2','de_halfmax_ell_maj']
    pix=206265./50000/.3 # from pc to arcsec and then to pix
    pix=pix*2.354 # new convention is to return the major axis for ell
    
    #k=['e_maj2','de_halfmax_ell_maj']
    pl.plot(t[k[0]][tfound]*pix,moments[k[1]][tindex]/(t[k[0]][tfound]*pix),'.') 
    pl.ylabel("new "+k[1]+" / old "+k[0])
    
    pl.xlabel("old "+k[0])
    
    delt=pl.absolute(t[k[0]][tfound]*pix-moments[k[1]][tindex])/(t[k[0]][tfound]*pix+moments[k[1]][tindex])*2
    z=pl.where(delt>0.1)[0]
    
    for i in range(len(z)):
        pl.text(t[k[0]][tfound][z[i]]*pix,moments[k[1]][tindex][z[i]]/(t[k[0]][tfound][z[i]]*pix),moments['id'][tindex][z[i]])
    #    pl.text(t[k[0]][z[i]]*pix,moments[k[1]][tindex][z[i]],moments['id'][tindex][z[i]])
    
    pl.savefig("pyprops."+datestr+"."+root+"old.new."+k[0]+".png")
        


if doplots:
#    pl.clf()
#    pl.plot(t['e_maj2']*pix,moments['de_halfmax_ell_maj'][tindex],'.')
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
    
    z=pl.where(pl.isnan(moments['de_halfmax_ell_maj']))[0]
    pl.plot(moments['halfmax_ell_maj'][z],moments['halfmax_ell_min'][z],'s',label='dec=nan')
    z=pl.where((moments['de_halfmax_ell_min']==0)*(moments['de_halfmax_ell_maj']>0))[0]
    pl.plot(moments['halfmax_ell_maj'][z],moments['halfmax_ell_min'][z],'cd',label='dec min=0')
    z=pl.where((moments['de_halfmax_ell_min']==0)*(moments['de_halfmax_ell_maj']==0))[0]
    pl.plot(moments['halfmax_ell_maj'][z],moments['halfmax_ell_min'][z],'r*',label='dec both=0')
    pl.legend(loc="best",prop={"size":10},numpoints=1)
    
    pl.savefig("pyprops."+datestr+"."+root+".measured_ellipses.png")
    
    
    pl.clf()
    pl.plot(moments['posang'],moments['halfmax_ell_min'],'.')
    pl.xlabel("posang")
    pl.ylabel("measured minor ell @halfmax")
    pl.plot([bm_pix[2],bm_pix[2]],pl.ylim(),'k',linestyle='dashed')
    pl.plot(pl.xlim(),[bm_pix[1],bm_pix[1]],'k',linestyle='dotted')
    pl.plot(pl.xlim(),[bm_pix[1]*1.12,bm_pix[1]*1.12],'k',linestyle='dotted')
    
    z=pl.where(pl.isnan(moments['de_halfmax_ell_maj']))[0]
    pl.plot(moments['posang'][z],moments['halfmax_ell_min'][z],'s',label='dec=nan')
    z=pl.where((moments['de_halfmax_ell_min']==0)*(moments['de_halfmax_ell_maj']>0))[0]
    pl.plot(moments['posang'][z],moments['halfmax_ell_min'][z],'cd',label='dec min=0')
    z=pl.where((moments['de_halfmax_ell_min']==0)*(moments['de_halfmax_ell_maj']==0))[0]
    pl.plot(moments['posang'][z],moments['halfmax_ell_min'][z],'r*',label='dec both=0')
    pl.legend(loc="best",prop={"size":10},numpoints=1)
    
    pl.savefig("pyprops."+datestr+"."+root+".measured_ellipses_angles.png")
    
    
    
    p  =moments['posang'].copy()
    dp=moments['de_posang'].copy()
    p=(p-bm_pix[2]+pl.pi/2)%pl.pi +bm_pix[2]-pl.pi/2
    dp=(dp-bm_pix[2]+pl.pi/2)%pl.pi +bm_pix[2]-pl.pi/2
    
    pl.clf()
    pl.plot(p,dp,'.')
    pl.plot([bm_pix[2],bm_pix[2]],pl.ylim(),'k')
    pl.plot(pl.xlim(),[bm_pix[2],bm_pix[2]],'k')
    pl.xlabel("meas posang")
    pl.ylabel("deconv posang")
    
    pl.plot(pl.xlim(),[bm_pix[2]+pl.pi/2,bm_pix[2]+pl.pi/2],'k',linestyle='dotted')

    pl.savefig("pyprops."+datestr+"."+root+".dec_ellipses_angles.png")
