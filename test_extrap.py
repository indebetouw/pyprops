import pickle
# pickle.dump([t,mom2x,mom2y],open("test_extrap.pkl","wb"))
t,mom2x,mom2y=pickle.load(open("test_extrap.pkl","rb"))

import pylab as pl
pl.ion()
pl.clf()

pl.plot(t,mom2x,'k')
tmax=t.max()

levs=[100,200,300,400,494]
c=['b','g','c','y','m']
for i in range(len(levs)):
    il=levs[i]
    e=extrap(t[0:il],mom2x[0:il])
    pl.plot([0,tmax],[e,0],color=c[i])
#    e2=extrap(t[0:il],mom2x[0:il],square=True)
    e2=extrap(t[0:il],mom2x[0:il],weight=False)
    pl.plot([0,tmax],[e2,0],color=c[i],linestyle='dashed')
    print e,e2
