# FU UVVIS
This is a collection of utilities for the Ocean Insight Spectrometers

## Prerequisites
```
pip install matplotlib numpy seabreeze
```

## Installation
```
git clone https://github.com/padalev/FU-UVVIS
cd FU-UVVIS
pip install .
```

## Usage
### Initialization
```
m = UVVIS.Measurement( integrationtime = 4000, averages = 30)
```

### Scope
```
m.scope()
```
Will show a live plot of the bare spectrum with averages applied.

### Set Dark
```
m.setdark()
```
Live plot of spectrum will open for reference. Confirm by clicking anywhere on the plot or press a key.

### Set Light
```
m.setlight()
```
Live plot of spectrum will open for reference. Confirm by clicking anywhere on the plot or press a key.

### Absorbance Measurement
```
m.absorbance(name)
```
Live plot of spectrum and absorbance will open for reference. Confirm by clicking anywhere on the plot or press a key.
Data will be saved to file beginning with the sample name.
Requires set dark and light spectra.

### Set Averages
```
m.setaverages(n)
```
Number of averages can be changed any time.
Dark and Light spectrum will have to be set again afterwards before actual measurement.

### Set Integration Time
```
m.setintegrationtime(t = None, auto = False)
```
Integration time can be set at any time in units of microseconds.
Dark and Light spectrum will have to be set again afterwards before actual measurement.
If `t` is set to a number it will be used as the new integration time.
If `auto` is set to `True` the integration time will be determined automatically by maximizing the spectral maximum.
Else a window will open for manual integration time setting (todo).

## Example
The following code will first ask for a light spectrum, then for a dark spectrum. Then three absorbance measurements of the samples are made using the list of sample names.

```
import UVVIS

m = UVVIS.Measurement(4000,30)
m.setlight()
m.setdark()
samples = ['sample_0','sample_1','sample_2']
for sample in samples:
  m.absorbance(sample)
```