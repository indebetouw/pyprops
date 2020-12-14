# what do you get if you calculate the second moment down to the half-max
# contour?
# for ellfit(halfmax pixels, wt=None),
# the variance is
# int_(-r)^r x^2 dx int_(-sqrt(r^2-x^2))^(sqrt(r^2-x^2)) dy / pi r^2
# = r^2/4.
# so ellfit returns a radius = r/2, where r is the extent of the half-intensity
# circle, r= sigma sqrt(2 ln(2))
# ellfit returns sigma sqrt(ln(2)/2) = 0.5887 sigma
#
# for ellfit(halfmax pixels, wt=I), its more complicated
# var = int_(-r)^r x^2 dx int_(-sqrt(r^2-x^2))^(sqrt(r^2-x^2)) exp(-0.5(x^2+y^2)) dy / pi r^2

from scipy.special import erf
import numpy as np

n=5000
uu=-1+2.*np.arange(n+1)/n # u=x/r r= sigma*sqrt(2*ln(2))
vv=-1+2.*np.arange(n+1)/n

# just y integral done, this is the integrand for the x integral

yint= erf(np.sqrt(np.log(2)*(1-uu**2))) - erf(-np.sqrt(np.log(2)*(1-uu**2)))

import matplotlib.pyplot as pl
pl.ion()
pl.clf()
pl.plot(uu,yint)

xyintegrand=np.sqrt(np.pi/2) *(2*np.log(2))**2 *uu**2 *np.exp(-np.log(2)*uu**2)
# times sigma^3

pl.clf()
du=2./n
pl.plot(uu,np.cumsum(xyintegrand)*du)

# normalization factor

xyintegrand0=np.sqrt(np.pi/2) *(2*np.log(2))**0.5 *np.exp(-np.log(2)*uu**2)
# times sigma

print(np.sum(xyintegrand)*du, np.sum(xyintegrand0)*du)
print(np.sqrt( np.sum(xyintegrand)/np.sum(xyintegrand0)) ) # times sigma

# get 0.671 sigma
