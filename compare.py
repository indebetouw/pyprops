def fixsizes(moments,bm_pix):
    # plausible code to impose lower limits on deconvolved ell sizes:

    # set minima on errors:
    z=pl.where( (pl.isnan(moments['de_halfmax_ell_maj']) | # failed
                 ((moments['de_halfmax_ell_maj']==0)*
                  (moments['de_halfmax_ell_min']==0)) | # both axes too small
                 (moments['de_halfmax_ell_maj']<(bm_pix[0]/2))) *
                (moments['dde_halfmax_ell_maj']<bm_pix[0]/2 ) )[0] 
    moments['dde_halfmax_ell_maj'][z]=bm_pix[0]/2
    z=pl.where( (pl.isnan(moments['de_halfmax_ell_maj']) | # failed
                 ((moments['de_halfmax_ell_maj']==0)*
                  (moments['de_halfmax_ell_min']==0)) | # both axes too small
                 (moments['de_halfmax_ell_maj']<(bm_pix[0]/2))) *
                (moments['dde_halfmax_ell_min']<bm_pix[1]/2 ) )[0] 
    moments['dde_halfmax_ell_min'][z]=bm_pix[1]/2

    # now fix the actual values
    z=pl.where( pl.isnan(moments['de_halfmax_ell_maj']) | # failed
                ((moments['de_halfmax_ell_maj']==0)*
                 (moments['de_halfmax_ell_min']==0)) | # both axes too small
                (moments['de_halfmax_ell_maj']<(bm_pix[0]/2)) )[0] # succeeded but < half beam
    # set to half beam:
    moments['de_halfmax_ell_maj'][z]=bm_pix[0]/2
    moments['de_halfmax_ell_min'][z]=bm_pix[1]/2
    moments['de_posang'][z]=bm_pix[2]


    # now deal with those that were only too small in the minor axis
    z=pl.where( (moments['de_halfmax_ell_min']<(bm_pix[1]/2))*
                (moments['dde_halfmax_ell_min']<(bm_pix[1]/2)) )[0] 
    moments['dde_halfmax_ell_min'][z]=bm_pix[1]/2
    z=pl.where( moments['de_halfmax_ell_min']<(bm_pix[1]/2) )[0] 
    moments['de_halfmax_ell_min'][z]=bm_pix[1]/2


import pickle
hcop,bm_pix,mhcop  =pickle.load(open("pyprops.20170105.dor.hcop.pkl",'rb'))
hcn,bm_pix,mchcn   =pickle.load(open("pyprops.20170105.dor.hcn.pkl",'rb'))
co13,bm_pix,mc13co =pickle.load(open("pyprops.20170106.dor.13co.pkl",'rb'))
co12,bm_pix,mc12co =pickle.load(open("pyprops.20170105.dor.12co.pkl",'rb'))
cs,bm_pix,mccs     =pickle.load(open("pyprops.20170105.dor.cs.pkl",'rb'))
# on edge:
z=[12,23,24,26,35,37,38,76,121]
for i in z:
    zz=pl.where(co13['id']==i)[0]
    co13['flux'][zz]=0.

fixsizes(hcn,bm_pix)
fixsizes(hcop,bm_pix)
fixsizes(co13,bm_pix)
fixsizes(cs,bm_pix)
fixsizes(co12,bm_pix)

pickle.dump([co12,co13,hcop,hcn,cs,bm_pix],open("pyprops.20170106.dor.fixed.pkl","wb"))

pl.subplots_adjust(left=.1,right=.95,bottom=.1,top=.95,wspace=.15,hspace=.15)
#=============================================================

pl.clf()
pl.errorbar(hcop['de_halfmax_ell_maj'],hcn['de_halfmax_ell_maj'],
            xerr=hcop['dde_halfmax_ell_maj'],yerr=hcn['dde_halfmax_ell_maj'],
            fmt='.')
pl.errorbar(hcop['de_halfmax_ell_min'],hcn['de_halfmax_ell_min'],
            xerr=hcop['dde_halfmax_ell_min'],yerr=hcn['dde_halfmax_ell_min'],fmt='.')
pl.xlabel("hco+ size")
pl.ylabel("hcn size")

pl.plot(pl.xlim(),pl.xlim(),'k')


#=============================================================
pl.clf()
hcopa=hcop['de_halfmax_ell_maj']*hcop['de_halfmax_ell_min']/bm_pix[0]/bm_pix[1]
dhcopa=hcopa*pl.sqrt( (hcop['dde_halfmax_ell_maj']/hcop['de_halfmax_ell_maj'])**2+
                     (hcop['dde_halfmax_ell_min']/hcop['de_halfmax_ell_min'])**2 )
hcna = hcn['de_halfmax_ell_maj']*hcn['de_halfmax_ell_min']/bm_pix[0]/bm_pix[1]
dhcna= hcna*pl.sqrt( (hcn['dde_halfmax_ell_maj']/hcn['de_halfmax_ell_maj'])**2+
                     (hcn['dde_halfmax_ell_min']/hcn['de_halfmax_ell_min'])**2 )
co13a = co13['de_halfmax_ell_maj']*co13['de_halfmax_ell_min']/bm_pix[0]/bm_pix[1]
dco13a= co13a*pl.sqrt( (co13['dde_halfmax_ell_maj']/co13['de_halfmax_ell_maj'])**2+
                     (co13['dde_halfmax_ell_min']/co13['de_halfmax_ell_min'])**2 )
z=pl.where( (co13['flux']>0)*(hcop['flux']>0) )[0]
pl.errorbar(hcopa[z],co13a[z],xerr=dhcopa[z],yerr=dco13a[z],fmt='.')
pl.xlabel("hco+ area")
pl.ylabel("13co area")

pl.plot(pl.xlim(),pl.xlim(),'k')
pl.xscale("log",nonposx='clip')
pl.yscale("log",nonposy='clip')
pl.xlim(.1,5)
pl.ylim(.1,5)
pl.plot(pl.xlim(),pl.xlim(),'k')

diff=pl.absolute(hcopa-co13a)*2./(hcopa+co13a)
z=pl.where( (diff>.8)*(co13a>co13a.min())*(co13['flux']>0)*(hcop['flux']>0) )[0]
for i in range(len(z)):
    pl.text(hcopa[z[i]],co13a[z[i]],hcop['id'][z[i]])


#14: 13CO really is larger, measurments probably ok.
#24: 13CO edge of field, both low SN
#27: 13CO really larger, less well resolved from neighbor clumps
#50: hco+ lower S/N, probably accurately measured as larger/blobbier
#85: 13CO low S/N, measurement probably not very good. (has large error bar)

pl.savefig("compare.hcop.13co.sizes.png")

#pl.clf()
#k='max'
#pl.errorbar(hcop[k],co13[k],xerr=hcop['d'+k],yerr=co13['d'+k],fmt='.')
#pl.xlabel("hco+ "+k)
#pl.ylabel("co13 "+k)


#=============================================================
# hcop already normalized by bmaj*bmin
hcopfact=pl.sqrt( hcopa/(1+hcopa) )
co13fact=pl.sqrt( co13a/(1+co13a) )
#pl.plot(hcop[k]/hcopfact,co13[k]/co13fact,'.')
pl.ylim(.01,10)

pl.clf()
k='max'
rat=hcop[k]/co13[k]
drat=rat*pl.sqrt( (hcop['d'+k]/hcop[k])**2 + (co13['d'+k]/co13[k])**2 )
de_rat=hcop[k]/co13[k]/hcopfact*co13fact
dde_rat=de_rat*pl.sqrt( (hcop['d'+k]/hcop[k])**2 + (co13['d'+k]/co13[k])**2 )

z=pl.where( (co13['flux']>0)*(hcop['flux']>0) )[0]
pl.errorbar(rat[z],de_rat[z],xerr=drat[z],yerr=dde_rat[z],fmt='o')
pl.xlabel("HCO+/13CO measured")
pl.ylabel("HCO+/13CO deconvolved")
pl.xscale("log",nonposx='clip')
pl.yscale("log",nonposy='clip')
pl.xlim(.02,.5)
pl.ylim(.02,.5)
pl.plot(pl.xlim(),pl.xlim(),'k')
pl.plot(pl.xlim(),pl.array(pl.xlim())*1.2,'k',linestyle='dotted')
pl.plot(pl.xlim(),pl.array(pl.xlim())*0.8,'k',linestyle='dotted')

pl.savefig("compare.hcop_13co.peakratio_deconv_20pct.png")




#=============================================================
hcnfact=pl.sqrt( hcna/(1+hcna) )

pl.clf()
rat=hcn[k]/hcop[k]
drat=rat*pl.sqrt( (hcn['d'+k]/hcn[k])**2 + (hcop['d'+k]/hcop[k])**2 )
de_rat=hcn[k]/hcop[k]/hcnfact*hcopfact
dde_rat=de_rat*pl.sqrt( (hcn['d'+k]/hcn[k])**2 + (hcop['d'+k]/hcop[k])**2 )

z=pl.where( (hcn['flux']>0)*(hcop['flux']>0) )[0]
pl.errorbar(rat,de_rat,xerr=drat,yerr=dde_rat,fmt='o')
pl.xlabel("HCN/HCO+ measured")
pl.ylabel("HCN/HCO+ deconvolved")
pl.xscale("log",nonposx='clip')
pl.yscale("log",nonposy='clip')
pl.xlim(.1,2)
pl.ylim(.1,2)
pl.plot(pl.xlim(),pl.xlim(),'k')
pl.plot(pl.xlim(),pl.array(pl.xlim())*1.2,'k',linestyle='dotted')
pl.plot(pl.xlim(),pl.array(pl.xlim())*0.8,'k',linestyle='dotted')

diff=pl.absolute(rat-de_rat)*2./(rat+de_rat)
z=pl.where( (diff>.5)*(co13['flux']>0)*(hcop['flux']>0) )[0]
for i in range(len(z)):
    pl.text(rat[z[i]],de_rat[z[i]],hcop['id'][z[i]])

pl.savefig("compare.hcn_hcop.peakratio_deconv_20pct.png")



#=============================================================
pl.clf()
#pl.plot(hcnfact,diff,'.')
#pl.xlabel("beam deconv. factor HCN")
pl.plot(hcn['de_halfmax_ell_maj']/hcn['dde_halfmax_ell_maj'],diff,'.')
pl.xlabel("HCN ellmax maj s.n")
pl.ylabel("HCN/HCO+ ratio, diff with/out deconv.")




pl.clf()
rat=hcna/hcopa
raterr=rat*pl.sqrt( (dhcna/hcna)**2+ (dhcopa/hcopa)**2 )
I=hcn['flux']
dI=hcn['dflux']
z=pl.where ((hcn['flux']>0)*(hcop['flux']>0) )[0]
pl.plot(I[z],rat[z],'x',label="HCN unresolved")
z=pl.where( (hcna>0.25)*(raterr<rat)*(hcn['flux']>0)*(hcop['flux']>0) )[0]
if len(z)>0:
    pl.errorbar(I[z],rat[z],xerr=dI[z],yerr=raterr[z],fmt='.b')

z=pl.where((hcna>0.25)*(hcn['flux']>0)*(hcop['flux']>0))[0]
pl.plot(I[z],rat[z],'o',label="HCN resolved")
pl.xscale("log",nonposx='clip')
pl.yscale("log",nonposy='clip')
pl.xlabel("HCN flux")
pl.ylabel("HCN size / HCO+ size")
pl.ylim(.05,10)
pl.legend(loc="best",prop={"size":10},numpoints=1)
z=pl.where( (rat>3) * (hcna>0.25) * (hcn['flux']>0) *(hcop['flux']>0) )[0]
for i in z:
    pl.text( I[i],rat[i],hcop['id'][i])
z=pl.where( (rat<.3333) * (hcna>0.25) * (hcn['flux']>0) *(hcop['flux']>0) )[0]
for i in z:
    pl.text( I[i],rat[i],hcop['id'][i])

pl.savefig("compare.hcn_hcop_sizes_flux.png")





pl.clf()
rat=hcopa/co13a
raterr=rat*pl.sqrt( (dco13a/co13a)**2+ (dhcopa/hcopa)**2 )
I=hcop['flux']
dI=hcop['dflux']
z=pl.where( (co13['flux']>0)*(hcop['flux']>0) )[0]
pl.plot(I[z],rat[z],'x',label="HCO+ unresolved")
z=pl.where( (hcopa>0.25)*(raterr<rat)*(co13['flux']>0) )[0]
if len(z)>0:
    pl.errorbar(I[z],rat[z],xerr=dI[z],yerr=raterr[z],fmt='.b')

z=pl.where((hcopa>0.25)*(co13['flux']>0)*(co13['flux']>0)*(hcop['flux']>0))[0]
pl.plot(I[z],rat[z],'o',label="HCO+ resolved")
pl.xscale("log",nonposx='clip')
pl.yscale("log",nonposy='clip')
pl.xlabel("HCO+ flux")
pl.ylabel("HCO+ size / 13CO size")
pl.xlim(1,300)
pl.ylim(.1,30)
pl.legend(loc="best",prop={"size":10},numpoints=1)
z=pl.where( (rat>3) * (hcopa>0.25) * (co13['flux']>0) *(hcop['flux']>0) )[0]
for i in z:
    pl.text( I[i],rat[i],hcop['id'][i])
z=pl.where( (rat<.3333) * (hcopa>0.25) * (co13['flux']>0) *(hcop['flux']>0) )[0]
for i in z:
    pl.text( I[i],rat[i],hcop['id'][i])

pl.savefig("compare.hcop_co13_sizes_flux.png")
