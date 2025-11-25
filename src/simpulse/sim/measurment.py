# sim/measurement.py

import numpy as np
import math as m


class MeasurementMixin:
    """
    Mixin holding all SNR / flux measurement methods.
    Injected into the Spectra class.
    """

    def write_snr(self):
        """Harry's fscrunch and L2 snr script"""
        base2 = self.burst_dedispersed
        quadsn = L2_clean(base2)
        fwhm = (m.sqrt(8.0 * m.log(2.0))) * self.width
        return f"{self.dm};{self.width};{fwhm};{quadsn}\n", quadsn

    def write_flux(self):
        """Compute L2_flux of the dedispersed burst."""
        base2 = self.burst_dedispersed
        flux = L2_flux(base2)
        return flux

def simulate(array, std=18, base=127, outtype=np.uint8):
    bkg = np.random.randn(array.shape[0], array.shape[1]) * std + base
    imprint = (bkg + array).astype(outtype)
    return imprint


def quick_snr(sf):
    # print("snr_check")
    ## single band snr method
    return np.sum(sf[sf > 0] ** 2) ** 0.5


def quad_sum(sf):
    # print("snr_check")
    ## multi-band snr method
    return np.sum(sf ** 2) ** 0.5


def fscrunch(array, prepost=2000):
    # this is an fscrunch which trims edges
    # directly copied from dynspec.py
    dedata = array.astype(np.float64)
    fscr = np.mean(dedata, axis=0)
    med = np.median(fscr)
    mad = np.median(np.abs(fscr - med))
    mask = (fscr - med) / mad < 5
    return fscr * mask


def L2_snr(base2):
    """Harry's fscrunch and L2 snr script"""
    simdata = simulate(base2, outtype=np.float64)  # base2 is the clean burst array
    fscrunched = np.sum((simdata.astype(np.float64)), axis=0)
    fscrun_mean = np.mean(fscrunched)
    fscrun_median = np.median(fscrunched)
    fscrun_mad = np.median(np.abs(fscrunched - fscrun_mean))  ##use MAD

    mask = np.sum(base2, axis=0) / fscrun_mad > 1  # find where pulse is after fscrunch
    sf = ((fscrunched - fscrun_median) / fscrun_mad)[mask]

    quadsn = (np.sum(sf ** 2) ** 0.5)
    return quadsn


def L2_clean(base2):
    """Harry's fscrunch and L2 snr script with no noise, assume rms/std is 1"""
    ydata = base2  # base2 is the clean burst array
    fscrunched = np.mean((ydata.astype(np.float64)), axis=0)
    mask = np.mean(base2, axis=0) > 0  # find where pulse is after fscrunch
    sf = fscrunched[mask]
    quadsn = (np.sum(sf ** 2) ** 0.5)
    return quadsn


def triangle_snr(base2):
    ## triangle method snr
    fscrunched = np.sum(base2, axis=0)
    slen = len(fscrunched)
    arr = []

    for i in range(slen):
        for j in range(i + 1, slen):
            tri = fscrunched[i:j]
            x = np.sum(np.abs(tri)) / np.sqrt(len(tri))
            arr.append(x)

    return np.max(arr)


def triangle_clean(base2):
    ## triangle method clean snr
    fscrunched = np.mean(base2, axis=0)
    slen = len(fscrunched)
    arr = []

    for i in range(slen):
        for j in range(i + 1, slen):
            tri = fscrunched[i:j]
            x = np.sum(np.abs(tri)) / np.sqrt(len(tri))
            arr.append(x)

    return np.max(arr)


def rollingbox(base2):
    """rolling boxcar filter"""
    fs = np.sum(base2, axis=0)
    l = len(fs)
    best = 0

    for w in range(1, l):
        box = np.convolve(fs, np.ones(w) / w, mode="same")
        sn = np.max(box)
        if sn > best:
            best = sn

    return best


def L2_flux(base2):
    """Harry's fscrunch and L2 snr script with no noise, assume rms/std is 1"""
    ydata = base2  # base2 is the clean burst array
    fscrunched = np.mean((ydata.astype(np.float64)), axis=0)
    mask = np.mean(base2, axis=0) > 0  # find where pulse is after fscrunch
    sf = fscrunched[mask]
    flux = np.sum(sf)
    return flux
