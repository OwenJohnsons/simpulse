"""
import_tests.py

Sanity checking imports after refactoring Harry's orginial repo. 
"""

print("=== IMPORT TESTS FOR SIMPULSE ===")

try:
    import simpulse
    print("PASS: import simpulse")
except Exception as e:
    print("FAIL: import simpulse -->", e)
    raise

try:
    from simpulse import Spectra, TimeSeries, fgrid
    print("PASS: from simpulse import Spectra, TimeSeries, fgrid")
except Exception as e:
    print("FAIL: public API imports -->", e)
    raise


try:
    from simpulse.sim.model import Spectra as S2
    from simpulse.sim.model import TimeSeries as T2
    from simpulse.sim.model import fgrid as F2
    print("PASS: simpulse.sim.model imports")
except Exception as e:
    print("FAIL: simpulse.sim.model imports -->", e)
    raise

try:
    import simpulse.sim.burst
    import simpulse.sim.noise
    import simpulse.sim.measurement
    print("PASS: simpulse.sim.* imports")
except Exception as e:
    print("FAIL: simpulse.sim.* imports -->", e)
    raise


try:
    from simpulse.io.fbio import makefilterbank
    print("PASS: simpulse.io.fbio imports")
except Exception as e:
    print("FAIL: simpulse.io.fbio -->", e)
    raise


print("\n=== FUNCTIONAL TEST ===")

try:
    m = Spectra()
    print("PASS: Spectra() instance created")
except Exception as e:
    print("FAIL: creating Spectra() -->", e)
    raise

try:
    m.create_filterbank("test_imports")
    print("PASS: create_filterbank()")
except Exception as e:
    print("FAIL: create_filterbank -->", e)
    raise

try:
    m.writenoise(100)
    print("PASS: writenoise()")
except Exception as e:
    print("FAIL: writenoise() -->", e)
    raise

try:
    orig, dd = m.burst(dm=100, width=1, A=10, nsamp=500)
    print("PASS: burst()  | shapes:", orig.shape, dd.shape)
except Exception as e:
    print("FAIL: burst() -->", e)
    raise

try:
    m.inject(dd)
    print("PASS: inject()")
except Exception as e:
    print("FAIL: inject() -->", e)
    raise


try:
    m.closefile()
    print("PASS: closefile()")
except Exception as e:
    print("FAIL: closefile() -->", e)
    raise


print("\n=== ALL TESTS PASSED ===")
