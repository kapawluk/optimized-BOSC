function ksstat=BOSC_compare_chi2(F,B,meanpower) 
% new function included with version 2.0
%
% This file is part of the Better OSCillation detection (BOSC) library.
%
% The BOSC library is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% The BOSC library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
%
% Copyright 2024 Kieran A. Pawluk, Tamari Shalamberidze and Jeremy B. Caplan.
%
% This function compares the power distribution at each frequency to the
% chi-square(2) distribution using KS tests
%
% parameters:
% F - vector containing frequencies sampled
% B - matrix containing power as a function of frequency (rows) and
%     time (columns). This is the time-frequency data.
% meanpower - geometric mean background power values by frequency aka
%             10^mean(log10(power)), estimated in newBOSC_bgfit.m
%
% returns:
% ksstat = KS statistic or KS d by frequency
%
B=sort(B,2); % sort power in each row (aka each frequency)
% scale power down to chi-square(2) with meanpower and the chi-square(2)
% geometric mean of e^(psi(1)+ln(2))
relpower=B./meanpower(:)*exp(psi(1)+log(2)); % scaled-down power values
% KS tests by frequency with scaled power and chi-square(2) distribution
[~,~,ksstat]=arrayfun(@(f) kstest(relpower(f,:),'CDF',makedist('Gamma','a',1,'b',2)),1:length(F));
% makedist.m is used to create a gamma(1,2) distribution, which is
% equivalent to a chi-square(2) distribution
