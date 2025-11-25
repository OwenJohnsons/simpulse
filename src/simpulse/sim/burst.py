# sim/burst.py

import numpy as np
import math as m
from scipy.signal import convolve

def dedisperse(dynamic_spectrum, dm, vif, fch1, tsamp):
    """Basic brute-force dedispersion."""
    out = np.zeros_like(dynamic_spectrum)
    nchan = dynamic_spectrum.shape[1]

    for i in range(nchan):
        delay = tidm(dm, vif[i], fch1)
        shift = int(delay / tsamp)
        out[:, i] = np.roll(dynamic_spectrum[:, i], -shift)

    return out

def boxcar_func(t, t0, a, width):
    y = np.zeros(t.shape[0])
    samp_diff = np.diff(t)[0]
    hw = width / 2
    p1 = np.argmin(np.abs(t - t0 + hw))
    y[p1:p1 + np.int64(width)] = a
    return y

class BurstMixin:
    """
    Mixin holding all burst-generation methods.
    This is injected into the Spectra class.
    """

    def burst(self,t0=100,dm=200,width=1,A=20,nsamp=5000,mode="boxcar",
              kscat=False,tau=0.1,alpha=4,offset=0.,dmoff=0,drift=0,bandfrac=None):
        """Create a dispersed pulse in noiseless data. Outputs both the dedispered and dedispersed pulse
        Parameters
        ----------
        mode : string
            Enter pulse shape used for injection: boxcar,scat,single
            boxcar: dynspec.boxcar_func
            scat: dynspec.spectra.scat_pulse_smear
            single: dynspec.spectra.single_pulse_smear
        width : float
            This is the 1-sigma of the gaussian, in units of ms.
            Note: for the boxcar it is the full width of the boxcar
        nsamp : int
            This sets the length of the array. Must be long enough for the dispersion track.
        A : float
            This is now the channel amplitude of the pulse with whichever mode, this parameter decides the injected value of the boxcar.
        """

        self.dm=dm
        self.width=width
        self.nsamp=nsamp
        tif=np.zeros((self.nchan*self.fbin,nsamp))
        tif2=np.zeros((self.nchan*self.fbin,nsamp))
        self.t0=t0

        if bandfrac is None:
            bandfrac = np.ones(self.nchan)

        ### time grid
        time = np.arange(nsamp) * self.tsamp

        ### compute frequency grid
        fgrid = self.vif.repeat(self.fbin)

        ### base arrays
        base = np.zeros((self.nchan*self.fbin, nsamp))
        ded = np.zeros((self.nchan*self.fbin, nsamp))

        ### injection loop
        for i in range(self.nchan*self.fbin):

            ### DM and drift delays
            tstart = (t0
                      + tidm(dm+dmoff, fgrid[i], self.fch1)
                      + pdrift(drift, fgrid[i], self.fch1)
                      + offset)

            ### scattering
            if kscat:
                tscat = tau * (fgrid[i]/1000)**(-alpha)

            ### choose shape
            if mode == "boxcar":
                pulse = boxcar_func(time, tstart, A, width)
            elif mode == "scat":
                pulse = scat_pulse_smear(time, tstart, width, A, tscat)
            elif mode == "single":
                pulse = single_pulse_smear(time, tstart, width, A)
            else:
                raise ValueError("Unknown mode {}".format(mode))

            base[i] = pulse

        ### Band fraction scaling
        ### reshape from (nchan*fbin, nsamp) to (nchan, fbin, nsamp)
        base = base.reshape(self.nchan, self.fbin, nsamp)
        base = base * bandfrac[:,None,None]
        base = base.reshape(self.nchan*self.fbin, nsamp)

        self.burst_original = base.reshape(self.nchan, self.fbin, nsamp).mean(1).T * bandfrac

        ### dedisperse
        self.burst_dedispersed = dedisperse(self.burst_original,
                                            dm=self.dm,
                                            vif=self.vif,
                                            fch1=self.fch1,
                                            tsamp=self.tsamp)
        return self.burst_original, self.burst_dedispersed

def single_pulse_smear(t, t0, width, A):
    """Single gaussian pulse"""
    return gaus_func(t, t0, width) * A

def gaus_func(t,t0,sigi):
    #ti=0### gaussian function returns a function with amplitude of 10
    A=1/np.sqrt(np.pi*2*(sigi**2))
    sit=A*np.exp(-1/2*((t-t0)**2)/(sigi**2)) ### model 0 in ravi 2018

    ### normalisation factor is 1/np.sqrt(np.pi*2*(sigi**2)) replace A with this term for a total of 1 pdf
    return sit


def scattering(t,t_0,tau1,alpha=4,v=1000):
    ###tau=tau1/1000 ## ms
    flux=np.zeros(len(t))
    flux[t>=t_0]=np.exp(-(t[t>=t_0]-t_0)/(tau1*(v/1000)**(-alpha)))
    return flux


def inverse_scattering(t,t_0,tau1,alpha=4,v=1000):
    ###tau=tau1/1000 ## ms
    flux=np.zeros(len(t))
    flux[t<t_0]=np.exp((t_0-t[t<t_0])/(tau1*(v/1000)**(-alpha)))
    return flux


def tidm(dm,vi,fch1):
    """ dispersion time delay offset """
    d=4148.808
    v=vi/1000
    top=fch1/1000
    dt = d * dm * (v**(-2) - top**(-2))
    return dt  ### ms


def pdrift(driftrate,vi,fch1):
    v=vi/1000
    top=fch1/1000
    dt = driftrate * (v - top)
    return dt
