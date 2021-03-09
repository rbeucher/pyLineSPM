import numpy as np

surface = np.zeros((100, 2))
surface[:, 0] = np.linspace(0, 8*np.pi, 100)
surface[:, 1] = np.cos(surface[:, 0])

precip_rate = np.ones(surface[:,0].shape)
erodibility = np.ones(surface[:,0].shape) * 0.1
uplift_rate = np.ones(surface[:,0].shape) * 0.01

def test_import():
    import pyLineSPM as spm

def test_spm_class():
    from pyLineSPM import WillettSPM
    spm = WillettSPM(surface, precip_rate, erodibility, uplift_rate)
    assert(len(spm.rivers) == 8)
    assert(np.allclose(surface, spm.surface))
    for river in spm.rivers:
        assert( 1.0 + river.surface[:,1].min() < 1e-2)
        assert( 1.0 - river.surface[:,1].max() < 1e-2)

    # Test get discharge
    discharge = spm.get_discharge()