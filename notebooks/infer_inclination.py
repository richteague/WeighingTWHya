"""
Script to infer the inclination of the disk assuming that the emission profile
is azimuthally symmetric.

Based on the Notebook by DFM.

Usage:

python infer_inclination.py /path/to/cube.fits
"""

import sys
import pymc3 as pm
import numpy as np
import exoplanet as xo
import theano.tensor as tt
from imgcube.cube import imagecube
np.random.seed(4321)

# Load up the data. Clip to reasonable bounds as we only want to fit out to 1".
# Note we convert to brightness temperature assuming the Rayleigh-Jeans law.
# This is just to make the numbers more manageable.

cube = imagecube(sys.argv[1], clip=1.1, kelvin='RJ')

# Pre-compute the meshgrid and flatten them.

X, Y = np.meshgrid(cube.xaxis, cube.yaxis)
X = X.flatten()
Y = Y.flatten()
y = np.array(cube.data.flatten(), dtype=np.float64)

# Take only positive values so we can work with the log and those within 1".

mask = (y > 0.0) & (np.hypot(X, Y) <= 1.0) & (np.random.rand(len(y)) <= 1.0)
X = X[mask]
Y = Y[mask]
y = np.log(y[mask])

# Define the deprojection functions.

def _deproject_coords(x, y, inc):
    return x, y / tt.cos(inc)

def _rotate_coords(x, y, PA):
    x_rot = x * tt.cos(PA) - y * tt.sin(PA)
    y_rot = y * tt.cos(PA) + x * tt.sin(PA)
    return x_rot, y_rot

def _get_cart_sky_coords(x0, y0):
    return X - x0, Y - y0

def _get_midplane_cart_coords(x0, y0, inc, PA):
    x_sky, y_sky = _get_cart_sky_coords(x0, y0)
    x_rot, y_rot = _rotate_coords(y_sky, x_sky, -PA)
    return _deproject_coords(x_rot, y_rot, inc)

def _get_midplane_polar_coords(x0, y0, inc, PA):
    x_mid, y_mid = _get_midplane_cart_coords(x0, y0, inc, PA)
    return tt.sqrt(y_mid**2 + x_mid**2), tt.arctan2(y_mid, x_mid)

# Make the mode and optimize the starting positions.

with pm.Model() as model:

    # Parameters
    x0 = pm.Normal("x0", mu=0.0, sd=0.1)
    y0 = pm.Normal("y0", mu=0.0, sd=0.1)
    inc = pm.Uniform("inc", lower=0.0, upper=0.5*np.pi)
    pa_deg = pm.Uniform("pa_deg", lower=0.0, upper=2.0*np.pi)
    pa = pa_deg * np.pi / 180.0
    mu = pm.Normal("mu", mu=np.mean(y), sd=100.0)
    slope = pm.Normal("slope", mu=0.0, sd=100.0)

    # Get the projected coordinates
    x_mid, y_mid = _get_midplane_cart_coords(x0, y0, inc, pa)
    r = tt.sqrt(y_mid**2 + x_mid**2)
    inds = tt.argsort(r)

    # Sort by radius
    r_sort = r[inds]
    y_sort = tt.as_tensor_variable(y)[inds]

    # Jitter & GP parameters
    logs2 = pm.Normal("logs2", mu=np.log(np.var(y)), sd=10)
    logw0_guess = np.log(2*np.pi/0.1)
    logtau = pm.Bound(pm.Normal, upper=0.0)(
        "logtau", mu=np.log(2*np.pi)-logw0_guess, sd=10)
    logw0 = np.log(2*np.pi) - logtau

    # We'll parameterize using the maximum power (S_0 * w_0^4) instead of
    # S_0 directly because this removes some of the degeneracies between
    # S_0 and omega_0
    logpower = pm.Normal("logpower",
                         mu=np.log(np.var(y))+4*logw0_guess,
                         sd=10)
    logS0 = pm.Deterministic("logS0", logpower - 4 * logw0)

    # Setup the GP
    kernel = xo.gp.terms.SHOTerm(log_S0=logS0, log_w0=logw0, Q=1/np.sqrt(2))
    gp = xo.gp.GP(kernel, r_sort, tt.exp(0.5*logs2) + np.zeros_like(y), J=2)

    # Compute the mean model
    line = mu + slope * r_sort

    # Compute the GP likelihood and predictions
    pm.Potential("loglike", gp.log_likelihood(y_sort - line))
    gp_pred = gp.predict() + line

    # Optimize to find the MAP
    map_soln = model.test_point
    map_soln = xo.optimize(map_soln, vars=[logs2, mu, slope])
    map_soln = xo.optimize(map_soln, vars=[x0, y0])
    map_soln = xo.optimize(map_soln, vars=[inc, pa_deg])
    map_soln = xo.optimize(map_soln, vars=[logs2, logpower, logtau])
    map_soln = xo.optimize(map_soln)

# Run the sampler.
np.random.seed(42)
sampler = xo.PyMC3Sampler(finish=200)
with model:
    sampler.tune(tune=2000, start=map_soln,
                 step_kwargs=dict(target_accept=0.9))
    trace = sampler.sample(draws=10000)

# Save the trace.
samples = pm.trace_to_dataframe(trace)
samples.to_pickle('%s' % (sys.argv[1].replace('.fits', '.trace.dat')))
