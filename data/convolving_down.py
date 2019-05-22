"""
Convolving down the images.
"""

# 230 GHz - 61 mas x 38 mas.
importfits(fitsimage='TWHya.230GHz.fits',
           imagename='TWHya.230GHz.im')
imsmooth(imagename='TWHya.230GHz.im',
         outfile='TWHya.230GHz.circ.im',
         kernel='gaussian',
         beam={'major': '0.062arcsec', 'minor': '0.062arcsec', 'pa': '0deg'},
         targetres=True)
exportfits(imagename='TWHya.230GHz.circ.im',
           fitsimage='TWHya.230GHz.circ.fits',
           dropstokes=True,
           velocity=True)
rmtables('TWHya.230GHz.im')
rmtables('TWHYa.230GHz.circ.im')

# 290 GHz - 37 mas x 26 mas.
importfits(fitsimage='TWHya.290GHz.fits',
           imagename='TWHya.290GHz.im')
imsmooth(imagename='TWHya.290GHz.im',
         outfile='TWHya.290GHz.circ.im',
         kernel='gaussian',
         beam={'major': '0.038arcsec', 'minor': '0.038arcsec', 'pa': '0deg'},
         targetres=True)
exportfits(imagename='TWHya.290GHz.circ.im',
           fitsimage='TWHya.290GHz.circ.fits',
           dropstokes=True,
           velocity=True)
rmtables('TWHya.290GHz.im')
rmtables('TWHYa.290GHz.circ.im')

# 345 GHz - 35 mas x 28 mas.
importfits(fitsimage='TWHya.345GHz.fits',
           imagename='TWHya.345GHz.im')
imsmooth(imagename='TWHya.345GHz.im',
         outfile='TWHya.345GHz.circ.im',
         kernel='gaussian',
         beam={'major': '0.036arcsec', 'minor': '0.036arcsec', 'pa': '0deg'},
         targetres=True)
exportfits(imagename='TWHya.345GHz.circ.im',
           fitsimage='TWHya.345GHz.circ.fits',
           dropstokes=True,
           velocity=True)
rmtables('TWHya.345GHz.im')
rmtables('TWHYa.345GHz.circ.im')

# 345 GHz - 35 mas x 28 mas.
importfits(fitsimage='TWHya.345GHz.fits',
           imagename='TWHya.345GHz.im')
imsmooth(imagename='TWHya.345GHz.im',
         outfile='TWHya.345GHz.largecirc.im',
         kernel='gaussian',
         beam={'major': '0.062arcsec', 'minor': '0.062arcsec', 'pa': '0deg'},
         targetres=True)
exportfits(imagename='TWHya.345GHz.largecirc.im',
           fitsimage='TWHya.345GHz.largecirc.fits',
           dropstokes=True,
           velocity=True)
rmtables('TWHya.345GHz.im')
rmtables('TWHYa.345GHz.largecirc.im')
