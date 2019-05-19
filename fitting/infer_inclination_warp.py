"""
Script to infer the inclination of the disk assuming that the emission profile
is azimuthally symmetric. Include a parameterised warp.

Based on the Notebook by DFM.

Usage:

python infer_inclination_warp.py /path/to/cube.fits
"""

import os
import sys
import pymc3 as pm
import numpy as np
import matplotlib.pyplot as plt

import exoplanet as xo
import theano.tensor as tt
from imgcube.cube import imagecube
np.random.seed(4321)

# Load up the data. Clip to reasonable bounds as we only want to fit out to 1".
# Note we convert to brightness temperature assuming the Rayleigh-Jeans law.
# This is just to make the numbers more manageable.

filepath = os.path.abspath(sys.argv[1])
filename = os.path.split(filepath)[1]
print(filepath, filename)
cube = imagecube(filepath, clip=1.1, kelvin='RJ')

if len(sys.argv) >= 3:
    nmax = int(sys.argv[2])
else:
    nmax = None

# Pre-compute the meshgrid and flatten them.

X, Y = np.meshgrid(cube.xaxis, cube.yaxis)
X = X.flatten()
Y = Y.flatten()
y = np.array(cube.data.flatten(), dtype=np.float64)

# Take only positive values so we can work with the log and those within 1".

mask = (y > 0.0) & (np.hypot(X, Y) <= 1.0)

if nmax is not None and mask.sum() > nmax:
    rands = np.random.rand(mask.sum())
    r = np.sort(rands)[nmax]
    mask[mask] = rands < r

X = X[mask]
Y = Y[mask]
y = np.log(y[mask])

# Define the deprojection function.

def deproject(x0, y0, inc, PA, w0, i0, r0):
    """
    Deproject the disk coordinates including a warped surface.
    If the disk is highly inclined or there are big perturbations,
    some pixels may be shadowed and convergence is not good.

    Args:
        x0 (float): Source center RA offset in arcsec.
        y0 (float): Source center Dec offset in arcsec.
        inc (float): Inclination of disk in radians.
        PA (float): Position angle of the disk in radians.
        w0 (float): Node angle of warp in radians.
        i0 (float): Magnitude of warp in radians.
        r0 (float): Radial scale of the warp in arcsec.

    Returns:
        rdisk (tensor array): Deprojected midplane radius in arcsec.
    """

    # Get the rotated disk coordiantes.
    x_sky, y_sky = X - x0, Y - y0
    x_rot = y_sky * tt.cos(PA) + x_sky * tt.sin(PA)
    y_rot = x_sky * tt.cos(PA) - y_sky * tt.sin(PA)

    # Iterate to account for warp. Increase iterations for better convergence.
    y_tmp = y_rot / tt.cos(inc)
    for _ in range(5):
        r_tmp = tt.sqrt(x_rot**2 + y_tmp**2)
        t_tmp = tt.arctan2(y_tmp, x_rot)
        w_tmp = i0 * tt.exp(-0.5*((r_tmp / r0)**2))
        z_tmp = r_tmp * np.tan(w_tmp) * np.sin(t_tmp - w0)
        y_tmp = y_rot / tt.cos(inc) - z_tmp * tt.tan(inc)
    return tt.sqrt(x_rot**2 + y_tmp**2)

# Make the mode and optimize the starting positions.

with pm.Model() as model:

    # Geometrical properties
    x0 = pm.Normal("x0", mu=0.0, sd=0.1)
    y0 = pm.Normal("y0", mu=0.0, sd=0.1)
    inc = pm.Uniform("inc", lower=0.0, upper=0.2*np.pi)
    pa_deg = pm.Uniform("pa_deg", lower=130., upper=170.)
    pa = pa_deg * np.pi / 180.0

    # Warp properties
    w0 = pm.Uniform("w0", lower=-0.5*np.pi, upper=0.5*np.pi)
    i0 = pm.Uniform("i0", lower=-0.3, upper=0.3)
    r0 = pm.Uniform("r0", lower=0.0, upper=3.0)

    # Mean model
    mu = pm.Normal("mu", mu=np.mean(y), sd=100.0)
    slope = pm.Normal("slope", mu=0.0, sd=100.0)

    # Get the projected coordinates
    r = deproject(x0, y0, inc, pa, w0, i0, r0)
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
    line = mu + slope * (r_sort - tt.mean(r_sort))

    # Compute the GP likelihood and predictions
    pm.Potential("loglike", gp.log_likelihood(y_sort - line))
    gp_pred = gp.predict() + line

    def optimize_geom(map_soln):
        map_soln = xo.optimize(map_soln, vars=[x0, y0])
        map_soln = xo.optimize(map_soln, vars=[inc_deg])
        map_soln = xo.optimize(map_soln, vars=[pa_deg])
        map_soln = xo.optimize(map_soln, vars=[inc_deg])
        map_soln = xo.optimize(map_soln, vars=[x0, y0])
        map_soln = xo.optimize(map_soln, vars=[inc_deg, pa_deg, x0, y0])
        map_soln = xo.optimize(map_soln, vars=[w0, i0, r0])
        map_soln = xo.optimize(map_soln, vars=[inc_deg, pa_deg, x0, y0,
                                               w0, i0, r0])
        return map_soln

    # Optimize to find the MAP
    map_soln = model.test_point
    map_soln = xo.optimize(map_soln, vars=[logs2, mu, slope])
    map_soln = optimize_geom(map_soln)
    map_soln = xo.optimize(map_soln, vars=[logs2, mu, slope])
    map_soln = optimize_geom(map_soln)
    map_soln = xo.optimize(map_soln, vars=[logs2, logpower, logtau])
    map_soln = optimize_geom(map_soln)
    map_soln = xo.optimize(map_soln)

    x_plot, y1_plot, y2_plot, y3_plot = xo.utils.eval_in_model(
        [r_sort, y_sort, line, gp_pred], map_soln)
    plt.plot(x_plot, y1_plot, ".k", label="data")
    plt.plot(x_plot, y2_plot, "r", label="linear model")
    plt.plot(x_plot, y3_plot, "g", label="gp model")
    plt.legend(fontsize=10)
    plt.title(filename)
    if nmax is None:
        fn = filename.replace('.fits', '.png')
    else:
        fn = filename.replace('.fits', '.{0}.png'.format(nmax))
    plt.savefig(fn, bbox_inches="tight")

print(map_soln)

# Run the sampler.
np.random.seed(42)
sampler = xo.PyMC3Sampler(finish=500)
with model:
    sampler.tune(tune=2000, start=map_soln,
                 step_kwargs=dict(target_accept=0.9))
    trace = sampler.sample(draws=10000)

# Save the trace.
samples = pm.trace_to_dataframe(trace)
if nmax is None:
    fn = filename.replace('.fits', '.warp.trace.dat')
else:
    fn = filename.replace('.fits', '.{0}.warp.trace.dat'.format(nmax))
samples.to_pickle(fn)
