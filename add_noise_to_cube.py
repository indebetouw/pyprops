def add_noise_to_cube(data, beamfwhm_pix, fluxmap=None):
    import pylab as pl
    pl.seed()
    s=data.shape
    noise=pl.randn(s[0],s[1],s[2])

    noisescale=1.
    if type(fluxmap)!=type(None):
        noisescale=1.26*fluxmap**2
        z=pl.where(pl.isnan(noisescale))
        if len(z[0])>0:
            noisescale[z]=1.
            

#    from astropy.convolution import convolve_fft,Gaussian2DKernel
#    psf=Gaussian2DKernel(stddev=beamfwhm_pix/2.354)
#    for i in range(s[0]):  # ASSUMES FIRST AXIS IS VEL
#        noise[i]=convolve_fft(noise[i]/noisescale,psf)#,interpolate_nan=True)

    from scipy.ndimage.filters import gaussian_filter
    for i in range(s[0]):  # ASSUMES FIRST AXIS IS VEL
        noise[i]=gaussian_filter(noise[i],beamfwhm_pix/2.354)/noisescale

        
    def mad(data, axis=None):
        return pl.nanmedian(pl.absolute(data - pl.nanmedian(data, axis)), axis)

    rms=mad(data) # rms of original cube
    current_rms=mad(noise)
    noise=rms*noise/current_rms # scale the noise to have the same rms as the data - there's a sqrt(2) problem I think
    
    return noise+data


