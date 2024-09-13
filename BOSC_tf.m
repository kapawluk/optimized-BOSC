function [B,T,F]=BOSC_tf(eegsignal,F,Fsample,wavenumber)
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
% Copyright 2010 Jeremy B. Caplan, Adam M. Hughes, Tara A. Whitten
% and Clayton T. Dickson.
%
% This function computes a continuous wavelet (Morlet) transform on
% a segment of EEG signal; this can be used to estimate the
% background spectrum (BOSC_bgfit) or to apply the BOSC method to
% detect oscillatory episodes in signal of interest (BOSC_detect).
%
% parameters:
% eegsignal - a row vector containing a segment of EEG signal to be
%             transformed
% F - a set of frequencies to sample (Hz)
% Fsample - sampling rate of the time-domain signal (Hz)
% wavenumber is the size of the wavelet (typically, width=6)
%	
% returns:
% B - time-frequency spectrogram: power as a function of frequency
%     (rows) and time (columns)
% T - vector of time values (based on sampling rate, Fsample) time of each
%     sample in seconds

st=1./(2*pi*(F/wavenumber));
A=1./sqrt(st*sqrt(pi));
B = zeros(length(F),length(eegsignal)); % initialize the time-frequency matrix
for f=1:length(F) % loop through sampled frequencies
  t=-3.6*st(f):(1/Fsample):3.6*st(f);
  m=A(f)*exp(-t.^2/(2*st(f)^2)).*exp(i*2*pi*F(f).*t); % Morlet wavelet
  y=conv(eegsignal,m); 
  y=abs(y).^2;
  B(f,:)=y(ceil(length(m)/2):length(y)-floor(length(m)/2));
end
T=(1:size(eegsignal,2))/Fsample;
