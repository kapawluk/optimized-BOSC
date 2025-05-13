    This file is part of the Better OSCillation detection (BOSC) library.

    The BOSC library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The BOSC library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

    Copyright 2010 Jeremy B. Caplan, Adam M. Hughes, Tara A. Whitten
    and Clayton T. Dickson.
    Updated 2024 Kieran A. Pawluk, Tamari Shalamberidze and Jeremy B. Caplan.


The optimized BOSC Method (BOSC version 2.0): Brief Tutorial
============================================================
* To use the recommended method of optimized BOSC, see the instructions below for the input variable "bgfitMethod" in the 	 
  BOSC_bgfit.m function.

Core functions and changes from version 1.0
-------------------------------------------
BOSC_tf.m - calculates time-frequency spectrograms based on the continuous Morlet wavelet transform
* No changes made

BOSC_bgfit.m - estimates the background spectrum
* Now includes different fitting options with the bgfitMethod parameter (string):
- bgfitMethod=“standard”: standard background fit in BOSC v1.0 with ordinary least-squares regression (polyfit.m)
- bgfitMethod=“robustfit”: exchanges polyfit.m for robustfit.m, a robust regression that downweighs outliers 
- bgfitMethod=“median”: BOSC v1.0 background fit is based on mean log10(power), this changes it to indirectly derive mean log10(power) from median log10(power), medians are less vulnerable to outliers
- bgfitMethod=“highpower”: performs initial polyfit.m regression, removes power values above threshold based on 99.9th percentile of chi-square(2) and results from 1st regression, then performs second polyfit.m regression without influence of high-power values
- bgfitMethod=“freqsubset”: performs initial polyfit.m regression, performs KS test on all frequencies, then performs second polyfit.m regression using frequencies with the lowest KS d values only (default: best 10 of 21)
- bgfitMethod=“optimized”: optimized BOSC v2.0 method, combines of “median”, “highpower” and “robustfit” options
	1) Initial bgfitMethod=“median” style regression
	2) Removes high-power values in the same manner as bgfitMethod=“highpower”
	3) Second regression again in bgfitMethod=“median” style but uses robustfit.m in place of polyfit.m
* Return variable meanpower is described as the estimated mean background power in BOSC v1.0 when it’s actually the estimated geometric mean (geomean) background power. Thresholds that scale power with a chi-square(2) using meanpower would previously also use the chi-square(2) mean=2 but now use the chi-square(2) geomean=e^(psi(1)+ln(2)).

BOSC_thresholds.m - calculates threshold values as a function of frequency, based on the estimate of the background (BOSC_bgfit). The power threshold is based on a percentile cutoff of the theoretical chi-square distribution of wavelet power values at each frequency. The duration threshold is converted into numbers of samples, scaling correctly for each frequency.
* Power threshold changed in the manner stated above since it uses meanpower, scaling done with chi-square(2) geomean=e^(psi(1)+ln(2)) instead of mean=2

BOSC_detect.m - detects oscillatory episodes in a target signal of interest, based on the thresholds calculated by BOSC_thresholds.m
* No changes made

BOSC_compare_chi2.m - compares power distribution at each sampled frequency with the chi-square(2) distribution using KS tests
* newly added function

Running the Oscillation Analyses
--------------------------------

example.m - a step-by-step walkthrough of how to run optimized BOSC method with annotations along the way, also contains code for various useful plots

simulate.m - creates simulated signals (to be optionally used in example.m) consisting of generated 1/f noise with oscillations in the form of sine waves added at certain times with varying frequency and amplitude

Important note: avoiding edge artifacts
---------------------------------------

As is important for data-analyses methods based on wavelet or other Fourier-domain method, data-windowing is an important element of time-frequency analysis. An appropriate data window is necessary to avoid edge artifacts due to not having sufficient sample points representing the signal when superimposing a particular wavelet used. For oscillation analysis, the duration threshold needs to be accounted for in addition to the wavelet window, both of which are frequency-dependent. This means a data window needs to be removed from both ends of the signal where the edges constrain the number of sample points up to the size of the window. A good rule of thumb is to select a data window, or “shoulder” that is (at a minimum) defined by the following: 

  (duration threshold + wavelet-window duration) X period of the slowest frequency sampled

In example.m, the necessary “shoulder” is expressed in numbers of samples:

	Fsample*(numcyclesthresh+width)/min(F)

For processes that don’t involve the duration threshold, such as the background fit, the data window or “shoulder” that needs to be stripped from both ends of the signal after the wavelet transform is:

  wavelet-window duration X period of the slowest frequency sampled

In example.m, this “background shoulder” is expressed in numbers of samples:

	Fsample*(width/min(F))
