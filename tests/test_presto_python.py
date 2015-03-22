import numpy as np
import newpresto as presto
import os

print "Testing FFT stuff...",
N = 20
x = np.random.standard_normal(N)
nx = presto.rfft(presto.rfft(x, -1), 1)
assert(np.allclose(x, nx, atol=1e-6))
print "success"

print "Testing FFTW call...",
cx = np.random.standard_normal(N).astype(np.complex64)
ncx = np.array(cx, copy=1)
presto.fftwcall(cx, -1)
presto.fftwcall(cx, 1)
assert(np.allclose(cx/N, ncx, atol=1e-6))
print "success"

print "Testing tablesixstepfft call...",
cx = np.random.standard_normal(N).astype(np.complex64)
ncx = np.array(cx, copy=1)
presto.tablesixstepfft(cx, -1)
presto.tablesixstepfft(cx, 1)
assert(np.allclose(cx/N, ncx, atol=1e-6))
print "success"

print "Testing reading infiles...",
x = presto.read_inffile("1937_DM71.02_zerodm.inf", verbose=False)
assert(x.telescope=="GBT")
assert(x.mjd_i==55267)
assert(x.dt==8.192e-05)
assert(x.numonoff==1)
assert(x.analyzer=="sransom")
print "success"

print "Testing writing infiles...",
x.analyzer="test"
x.name="xxx"
x.dt=0.125
presto.write_inffile(x, verbose=False)
y = presto.read_inffile("xxx", verbose=False)
assert(y.analyzer=="test")
assert(y.bary==0)
assert(y.numonoff==1)
assert(y.dt==0.125)
os.remove("xxx.inf")
print "success"

print "Testing allocation and freeing of memory...",
for ii in range(1024):
    a = presto.gen_fvect(1024 * 32768)
    del a
for ii in range(1024):
    a = presto.gen_cvect(1024 * 16384)
    del a
print "success"

print "Testing psrparams and orbitparams stuff...",
psr = presto.psrepoch("J0737-3039A", 56000.0, verbose=False)
assert(round(psr.dm-48.92, 7)==0)
assert(round(psr.orb.p-8834.534998272, 7)==0)
print "success"

print "Testing spectralpower and spectralphase...",
a = np.arange(5.0) + complex(0.0, 1.0)
assert(np.allclose(presto.spectralpower(a),
                   np.arange(5.0)**2.0 + 1))
assert(np.allclose(presto.spectralphase(a),
                   np.array([90., 45., 26.56505203, 18.43494797, 14.03624344])))
print "success"

print "Testing vector shifting / rotation...",
a = np.arange(4, dtype=np.float32)
presto.frotate(a, 1)
assert(np.allclose(a, np.array([1, 2, 3, 0])))
a = np.arange(4, dtype=np.float64)
presto.drotate(a, 1)
assert(np.allclose(a, np.array([1, 2, 3, 0])))
print "success"

print "Testing orbit integration stuff...",
orb = presto.orbitparams()
orb.p = 10000.0
orb.e = 0.1
orb.x = 1.0
orb.t = 1234.0
orb.w = 75.0
orb.pd = orb.wd = 0.0
E0 = presto.keplers_eqn(orb.t+0.0, orb.p, orb.e, 1e-15)
E1 = presto.keplers_eqn(orb.t+100.0, orb.p, orb.e, 1e-15)
E2 = presto.keplers_eqn(orb.t+200.0, orb.p, orb.e, 1e-15)
E3 = presto.keplers_eqn(orb.t+300.0, orb.p, orb.e, 1e-15)
Es = np.asarray([E0, E1, E2, E3])
Es_check = np.asarray([ 0.85050653, 0.9175909,
                        0.9842971, 1.05061346])
assert(np.allclose(Es, Es_check))
Es_new = presto.dorbint(E0, 4, 100.0, orb)
assert(np.allclose(Es_new, Es_check))
presto.E_to_v(Es, orb)
Vs_check = np.asarray([-112.15558594, -122.45299212,
                       -131.9991447, -140.76659065])
assert(np.allclose(Es, Vs_check))
minv, maxv = presto.binary_velocity(300.0, orb)
minv *= presto.SOL/1000.0
maxv *= presto.SOL/1000.0
assert(round(minv-Vs_check.min(), 7)==0)
assert(round(maxv-Vs_check.max(), 7)==0)
print "success"