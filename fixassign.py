# 79 is good in ch <=27, in ch28 we get 77,66, 78 to S
#       67   72    79   77   66   78      62
# ch28: ok   ok   big  not   w/78 w/66   faint
#   29  fnt tooS  ok    ok   w/78  ok    -  
# 77,79 poorly resolved in 13CO, HCN, ok in HCO+
# 66,78 also poorly resolved
# 118 could be larger
dir="/lustre/cv/users/rindebet/30dor/co.analysis/newcprops/"
f="hcop.20170101.clfind.assign.fits"
d=  fits.getdata(dir+f)
h=fits.getheader(dir+f)

d[36:38,220,238]=0
z=pl.where(d==31)
d[z]=32
z=pl.where(d==40)
d[z]=39

d[28,241,205:214]=72
d[29,237,211:217]=79
d[29,236,213:218]=79
d[29,235,215:219]=79
d[29,224:237,191:198]=0
d[29,219,202:207]=77
d[29,218,204:209]=77
d[29,217,206:213]=77
d[29,216,208:213]=77
d[29,215,210:214]=77
d[30,241,203:216]=72
d[30,236:239,212:217]=72
d[31,215:218,205:214]=77
d[31,226,206:208]=79
d[31,225:227,208:210]=79
d[32,240,201:204]=72
d[32,239,203:207]=72
d[32,238,206:210]=72
d[32,237,209:213]=72
d[32,236,212:216]=72
d[33,239,202:212]=72
d[33,238,204:215]=72
d[33,237,208:216]=72
d[33,236,211:216]=72
d[35,219:226,193:200]=77
d[35,219:226,192]=0
d[35,220,222:226]=0
d[35,227,201]=77
d[36,227,195:200]=77
d[36,226,195:202]=77
d[36,225,196:202]=77
d[36,224,198:202]=77
d[36,223,200:202]=77
d[36,222,200:202]=77


#27 is ok - it really is smaller in HCO+ than 13CO
#50 ok - ass could be larger - lower SN in HCO but maybe really larger
#93 ok - runs into 58 to the S, which is not well deliniated. 
# if we want 58 separately, it should be better assigned, and include 
# parts of 93.  so for now
z=pl.where(d==58)
d[z]=0

# 113 ass should be larger for hco+, but its not really htere in the other
# lines so it really is larger in hco+

# 24 - poorly separated from 23
z=pl.where(d==24)
d[z]=23

# 36 ok just not very bright


# now, split 68 and create 122
z=pl.where( (d[25]==52)|(d[25]==68) )
d[25,z[0],z[1]]=122
z=pl.where( (d[26]==68) )
d[26,z[0],z[1]]=122
d[26,229,255:258]=122
d[26,229:233,256]=122
d[26,231:235,257]=122
d[26,232:236,258:260]=122


z=pl.where( (d[27]==68) )
zz=pl.where(z[1]>248)[0]
d[27,z[0][zz],z[1][zz]]=122

z=pl.where( (d[28]==68) )
zz=pl.where(z[1]>249)[0] # x>249
d[28,z[0][zz],z[1][zz]]=122

z=pl.where( (d[29]==68) )
zz=pl.where((z[1]+z[0])>479)[0] # x+y>479
d[29,z[0][zz],z[1][zz]]=122

z=pl.where( (d[30]==68) )
zz=pl.where(z[1]>256)[0] 
d[30,z[0][zz],z[1][zz]]=0

z=pl.where( (d[32]==68) )
d[26,z[0],z[1]]=0


d[20:24,224:233,251:259]=52
d[35:38,233:237,175:180]=105


# go ahead and merge central clumps
z=pl.where(d==77)
d[z]=79
z=pl.where(d==78)
d[z]=66

fits.writeto(dir+"hcop.20170105.clfind.assign.fits",d,header=h,clobber=True)
