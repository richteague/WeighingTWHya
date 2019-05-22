# `/data/` Data Description

The data have been originally presented in two papers. The original Band 7 data was presented in [Andrews et al. (2016)](https://ui.adsabs.harvard.edu/#abs/2016ApJ...820L..40A/abstract) with images available from [Sean Andrew's website](https://www.cfa.harvard.edu/~sandrews/). [Huang et al (2018)](https://ui.adsabs.harvard.edu/#abs/2018ApJ...852..122H/abstract) combined several datasets to image TW Hya at both 230GHz (1.3mm) and 345GHz (870um), Band 6 and 7, respectively, as well as a combined dataset at 290GHz (1mm).

For consistency we take the three images from Huang et al. (2018) (however testing shows negligible difference from those presented in Andrews et al. 2016). As the images have non-circular beams, we also use `imsmooth` in CASA to convolve the image down to a circular beam.

A 12CO J = 2 - 1 image combining several projects worth of data was also presented in Huang et al. (2018). The data are available from the [Harvard Dataverse](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/PXDKBC). The 13CO J = 2 - 1, shown in the Appendix, is available upon request.

From these line images, rotation maps were generated using [`bettermoments`](https://github.com/richteague/bettermoments) with the standard values. For example:

```bash
$ bettermoments dataverse_TWHya_CO_cube.fits
```

which created both maps of the line center, `dataverse_TWHya_CO_cube.v0.fits`, and the uncertainties on the line center, `dataverse_TWHya_CO_cube.dv0.fits`.

### Continuum Data - `./data/cont/`

* `TWHya.230GHz.fits` - Band 6 data with a resolution of 60.7 mas x 37.9 mas.

* `TWHya.290GHz.fits` - combined Band 6 and Band 7 data with a resolution of 37.2 mas x 25.9 mas.

* `TWHya.345GHz.fits` - Band 7 data with a resolution of 35.1 mas x 27.5 mas.

* `TWHya.230GHz.circ.fits` - 230GHz image smoothed down to a circular beam of 62 mas.

* `TWHya.290GHz.circ.fits` - 290GHz image smoothed down to a circular beam of 38 mas.

* `TWHya.345GHz.circ.fits` - 345GHz image smoothed down to a circular beam of 36 mas.

* `TWHya.345GHz.largecirc.fits` - 345GHz image smoothed down to a circular beam of 62 mas to match the Band 6 data.

### Line Data - `./data/line/`

* `TWHya.12CO_v0.fits` - line center map of the 12CO data with a resolution of 139 mas x 131 mas.

* `TWHya.12CO_dv0.fits` - uncertainty on the line center of the 12CO data with a resolution of 139 mas x 131 mas.

* `TWHya.13CO_v0.fits` - line center map of the 13CO data with a resolution of 0.45" x 0.35".

* `TWHya.13CO_dv0.fits` - uncertainty on the line center of the 13CO data with a resolution of 0.45" x 0.35".

* `TWHya.12CO_M1.fits` - first moment map clipping at 5 sigma.

* `TWHya.12CO_dM1.fits` - uncertainty on the first moment map.
