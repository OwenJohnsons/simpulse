# sim/model.py

import numpy as np
import math as m
from simpulse.io.fbio import makefilterbank
from scipy.signal import convolve
import matplotlib.pyplot as plt

# Import mixins (implemented in other files)
from .noise import NoiseMixin
from .burst import BurstMixin
from .measurement import MeasurementMixin

def freq_splitter_idx(n, skip, end, bwchan, fch1):
    ### generates the frequency of channels and then group them into subbands, 
    ### also returns an array that records the channel numbers of each subband
    dw = (end - skip) / n
    vi = (np.arange(n) + 0.5) * dw * bwchan
    base = fch1 + skip * bwchan
    vi = base + vi
    chan_idx = (np.arange(n) * dw + skip)
    chan_idx = np.append(chan_idx, end).astype(np.int64)
    return vi, chan_idx


class TimeSeries:
    def __init__(self, tsamp=1, nsamp=1000, bins=10):
        """initiate function for creating a mock time series. This sets up the frequency.
        Parameters
        ----------
        tsamp : float
            time resolution (ms)
        nsamp : int
            This sets the length of the array. Must be long enough for scattering tail and dispersion track
        bins : int
            grid resolution of the array,
        """
        # self.fch=fch
        # self.bwchan=bwchan
        self.tsamp = tsamp
        self.nsamp = nsamp
        time = np.arange(nsamp) * tsamp
        matrix = np.ones((nsamp, bins)) * np.linspace(-0.5, 0.5, bins) * tsamp
        timematrix = (np.ones((nsamp, bins)).T * time).T
        finergrid = (matrix + timematrix).flatten()
        self.grid = finergrid
        self.x_time = time

    def boxcar(self, t0, width, a):
        tims = boxcar_func(self.x_time, t0, A, width)
        self.spectra = tims / np.max(tims) * a
        return self.spectra

    def pulse(self, t0, width, a):
        tims = np.mean(single_pulse(self.grid, t0, width, 100).reshape(self.nsamp, -1), axis=1)
        self.spectra = tims / np.max(tims) * a
        return self.spectra

    def scatp(self, t0, width, a, tau):
        tims = np.mean(scat_pulse(self.grid, t0, tau, width, 0, 100, 1000).reshape(self.nsamp, -1), axis=1)
        self.spectra = tims / np.max(tims) * a
        return self.spectra

    def inverse_scatp(self, t0, width, a, tau):
        tims = np.mean(invscat_pulse(self.grid, t0, tau, width, 0, 100, 1000).reshape(self.nsamp, -1), axis=1)
        self.spectra = tims / np.max(tims) * a
        return self.spectra



class Spectra(NoiseMixin, BurstMixin, MeasurementMixin):
    def __init__(self, fch1=1100, nchan=336, bwchan=1, tsamp=1,
                 nbits=8, fbin=10, tbin=10):
        """initiate function for creating a mock dynamic spectrum data. This sets up the header.
        Parameters
        ----------
        fch1 : float
            First channel centre frequency (MHz)
        nchan : int
            Total number of channels
        bwchan : float
            channel bandwidth (MHz)
        tsamp : float
            time resolution (ms)
        """

        self.fch1 = fch1
        self.nchan = nchan
        self.bwchan = bwchan
        self.tsamp = tsamp
        self.nbits = nbits
        self.fbin = fbin
        self.tbin = tbin

        # Frequency grid
        vi, chan_idx = freq_splitter_idx(nchan, 0, nchan, bwchan, fch1)
        self.vif = vi
        self.chan_idx = chan_idx

        # Header for filterbank writing
        self.header = {
            "telescope_id": 6,
            "fch1": fch1,
            "foff": -bwchan,
            "nchans": nchan,
            "tsamp": tsamp / 1000,
            "nbits": nbits,
        }

    def create_filterbank(self, file_name, std=np.sqrt(336), base=127):
        """Create a mock dynamic spectrum filterbank file.
        Parameters
        ----------
        file_name : string
            filename of filterbank
        std : float
            standard deviation of white noise, for normalised noise after fscrunching, set to sqrt(nchan)
        base : float
            base level of array
        """
        self.filterbank = makefilterbank(file_name + ".fil", header=self.header)
        self.fil_std = std
        self.fil_base = base

    def closefile(self):
        """Close writing filterbank"""
        self.filterbank.closefile()

    def inject(self, array):
        """Create a mock dynamic spectrum filterbank file.
        Parameters
        ----------
        array : numpy array object
            the burst array data to be injected into the filterbank object
        """
        scaledarray = array * self.fil_std / np.sqrt(array.shape[0])
        bkg = (np.random.randn(array.shape[0], array.shape[1]) *
               self.fil_std + self.fil_base)
        imprint = (bkg + scaledarray).astype(np.uint8)
        self.filterbank.writeblock(imprint)
        self.injected_array = imprint



class fgrid:
    def __init__(self, fch1=1000, bwchan=1, nchan=336, tsamp=1,
                 nsamp=1000, tbin=10, fbin=10):
        """Simulate a burst in a higher resolution grid. tgrid is the higher resolution 
        while fgrid is the final dynamic higher resolution
        Parameters
        ----------
        fch1 : float
            Upper frequency of array, not first channel frequency (MHz)
        bwchan: float
            Final product Channel bandwidth (MHz)
        nchan : int
            Number of channels of final product
        nsamp : int
            This sets the length of the array.
        tsamp : float
            time resolution (ms)
        tbin : int
            grid time resolution
        fbin : int
            grid frequency resolution
        """
        self.fch1 = fch1
        self.bwchan = bwchan
        self.nchan = nchan
        self.tsamp = tsamp
        self.nsamp = nsamp
        time = np.arange(nsamp) * tsamp

        tims = TimeSeries(tsamp=tsamp, nsamp=nsamp, bins=tbin)
        self.tims = tims
        self.tgrid = tims.grid
        self.x_time = tims.x_time

        vif, chan_idx = freq_splitter_idx(nchan, 0, nchan, bwchan, fch1)
        self.vif = vif
        self.chan_idx = chan_idx
        self.fbin = fbin
        self.tbin = tbin

        vif2, chan_idx2 = freq_splitter_idx(nchan * fbin, 0, nchan * fbin,
                                            bwchan / fbin, fch1 - bwchan * 0.5)
        self.fgrid = vif2

        self.array = np.zeros((nchan * fbin, nsamp))

    def pulse(self, t0, width, A, tau=10, alpha=4, dm=0, mode='gaussian', drift=0, dmerr=0):
        """Simulates pulse in datagrid
        Parameters
        ----------
        mode : string
            gaussian, scat, scat_r
        width : float
            Gaussian sigma (ms)
        dm : float
            dedispersed DM (no arrival time delay), for smearing simulation
        A : float
            amplitude
        """
        fbin = self.fbin
        for i in range(self.nchan):
            for j in range(fbin):
                tstart = (t0 + tidm(dm, self.fgrid[i*fbin+j], self.fgrid[i*fbin]) +
                          tidm(dmerr, self.fgrid[i*fbin+j], self.fch1) +
                          pdrift(drift, self.fgrid[i*fbin+j], self.fch1))

                smear = delta_t(dm + dmerr, self.fgrid[i*fbin+j], self.bwchan / fbin)
                smeared = np.sqrt(smear**2 + width**2)

                if mode == 'gaussian':
                    self.array[i*fbin+j] = self.tims.pulse(tstart, width, A)
                elif mode == 'scat':
                    tscat = tau * (self.fgrid[i*fbin+j] / 1000)**(-alpha)
                    self.array[i*fbin+j] = self.tims.scatp(tstart, width, A, tscat)

        self.model_burst = np.mean(self.array.reshape(self.nchan, fbin, self.nsamp), axis=1)
        return self.model_burst
    
