dir='/Users/remy/lustre/naasc/users/rindebet/30dor/chevance/'
datafile=dir+'dor12coS_r1.5.fapex.trim.fits'
fluxfile=dir+'dor12coS_r1.5.pb.fapex.ch0.trim.fits'
assignfile=dir+'20130410.assign_out.fits'
root="dor.rimc.12co"

import sys
gitpaths=['/Users/remy/lustre/naasc/users/rindebet/github/pyprops/']
for gitpath in gitpaths:
    if not gitpath in sys.path:
        sys.path.insert(0,gitpath)


from pyprops import pyprops

pyprops(datafile,fluxfile,assignfile,root,assignfile2=None,montecarlo=0,doplots=False):
