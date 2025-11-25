# sim/noise.py

import numpy as np


class NoiseMixin:
    """
    Noise-related methods mixed into Spectra.

    These methods assume that:
    - self.filterbank, self.fil_std, self.fil_base are set by create_filterbank()
    """

    def writenoise(self, nsamp=5000):
        """Write a block of white noise into the filterbank.
        Parameters
        ----------
        nsamp : int
            length of noise in units of tsamp
        """
        self.filterbank.writenoise(nsamp, self.fil_std, self.fil_base)
