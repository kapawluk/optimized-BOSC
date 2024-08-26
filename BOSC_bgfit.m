%    This file is part of the Better OSCillation detection (BOSC) library.
%
%    The BOSC library is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    The BOSC library is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
%
%    Copyright 2010 Jeremy B. Caplan, Adam M. Hughes, Tara A. Whitten
%    and Clayton T. Dickson.

function [pv,meanpower]=BOSC_bgfit(F,B,bgfitMethod)
% version 2.0
% This function estimates the background power spectrum via a regression 
% fit to the power spectrum in log-log coordinates using a selection of:
% normal BOSC, individual modifications or optimized BOSC
% 
% parameters:
% F - vector containing frequencies sampled
% B - matrix containing power as a function of frequency (rows) and
%     time (columns). This is the time-frequency data.
%
% returns:
% pv = contains the slope and y-intercept of regression line
% meanpower - geometric mean background power values by frequency aka
%             10^mean(log10(power)), estimated from regression
%
if bgfitMethod=="standard"
%%% normal BOSC %%%
pv=polyfit(log10(F),mean(log10(B),2),1); % linear regression
meanpower=10.^(polyval(pv,log10(F))); % transform back to natural units (power; usually uV^2/Hz)
elseif bgfitMethod=="robustfit"
%%% robustfit %%%
pv=flip(robustfit(log10(F),mean(log10(B),2)')); % robust regression
meanpower=10.^(polyval(pv,log10(F))); % transform back to natural units
elseif bgfitMethod=="median"
%%% median-based regression, log-means calculated from log-medians %%%
pv=polyfit(log10(F),(psi(1)+log(2)+median(log(B),2)-log(chi2inv(0.5,2)))/log(10),1); % linear regression
meanpower=10.^(polyval(pv,log10(F))); % transform back to natural units % transform back to natural units
elseif bgfitMethod=="highpower"
%%% remove high power %%%
    pv=polyfit(log10(F),mean(log10(B),2),1); % linear regression
    meanpower=10.^(polyval(pv,log10(F))); % transform back to natural units
    % threshold of 99.9th percentile of chi-square(2) (chi2inv(0.999,2)) * meanpower
    B(B>=(chi2inv(0.999,2)*meanpower(:)/exp(psi(1)+log(2))))=NaN; % power exceeding frequency-specific thresholds -> NaN
    pv=polyfit(log10(F),mean(log10(B),2,'omitnan'),1); % linear regression, NaNs excluded
    meanpower=10.^(polyval(pv,log10(F))); % transform back to natural units
elseif bgfitMethod=="freqsubset"
    %%% frequency subset with ks d %%%
    pv=polyfit(log10(F),mean(log10(B),2),1); % linear regression
    meanpower=10.^(polyval(pv,log10(F))); % transform back to natural units
    ks_d=BOSC_compare_chi2(F,B,meanpower); % perform KS tests to get KS d values
    [~,idx]=mink(ks_d,10); % use frequencies with lowest KS d (best 10, can be changed)
    pv=polyfit(log10(F(idx)),mean(log10(B(idx,:)),2),1); % linear regression using frequency subset
    meanpower=10.^(polyval(pv,log10(F))); % transform back to natural units
elseif bgfitMethod=="optimized"
    %%% optimized BOSC %%%
    % combination of median-based regression, high-power removal and robustfit
    pv=polyfit(log10(F),(psi(1)+log(2)+median(log(B),2)-log(chi2inv(0.5,2)))/log(10),1); % linear regression, log-medians -> log-means
    meanpower=10.^(polyval(pv,log10(F))); % transform back to natural units
    % threshold of 99.9th percentile of chi-square(2) (chi2inv(0.999,2)) * meanpower
    B(B>=(chi2inv(0.999,2)*meanpower(:)/exp(psi(1)+log(2))))=NaN; % power exceeding frequency-specific thresholds -> NaN
    pv=flip(robustfit(log10(F),(psi(1)+log(2)+median(log(B),2,'omitnan')-log(chi2inv(0.5,2)))/log(10)')); % robust regression, log-medians -> log-means
    meanpower=10.^(polyval(pv,log10(F))); % transform back to natural units
else
    fprintf("error: bgfitMethod invalid input\n");
end