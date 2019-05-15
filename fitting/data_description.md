# `/fitting/` Data Description

This folder, `/fitting/`, contains all the samples from the MCMC runs with various parameters.

## Accessing the Data

These can easily be read using `pandas`,

```python
import pandas as pd
samples = pd.read_pickle('path/to/trace.dat')
```

where `samples` is a dictionary containing all the posterior samples. For example,

```python
print('{:3f}'.format(np.median(samples['x0'])))
```

will print the median value of `x0`.

## File Descriptions

In general the files have the same format:

```
[TYPE].[FREQUENCY]GHZ.[circ.][DESCRIPTION.]trace.dat
```

where `[TYPE]` is the sort of data used, either `TWHya` for observations or `Model` for one of the model files in `../models/`, `[FREQUENCY]` is the frequency of the continuum used, either 230 GHz or 345 GHz. `[circ.]` indicates whether the continuum has been smoothed to a circular beam or not. Finally `[DESCRIPTION]` provides a cryptic description of the fitting:

| Desciption | Meaning |
| :-: | - |
| `[none]` | Default image with all values inwards of 1 arcsecond. |
| `inner0.7` | Pixels within the inner 0.7 arcseconds. This is roughly the 10sigma cut off. |
| `PAprior` | With a Gaussian prior provided for the position angle centered on 151 degrees with a standard deviation of 2 degrees. |
| `tightPAprior` | As above but with a standard deviation of 0.1 degrees. |
