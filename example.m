width=6; %  - the wavenumber of the Morlet wavelet (Grossman & Morlet, 1985)

Fsample=500; % Sampling rate of signal, in Hz.

F=(2^(1/4)).^(0:20); % Frequency sampling (resolution) for spectral analysis.  
% Logarithmic sampling of frequencies is recommended for computational purposes when using wavelets. The range can be modified but it must be within the bandpass of the amplifier and/or filter(s).

percentilethresh=.95; % Confidence level (range 0 to 1; percentile of the CDF divided by 100) of the estimated chi-square distribution of spectral power used for the power threshold. Increase to make more conservative.

numcyclesthresh=3; % Duration threshold expressed in numbers of oscillation cycles. Increase to make more conservative.

% empirical signal 
%bgsignal=EEG_PZ;

% simulated signal 
addpath(genpath('simulate'));
bgsignal=simulate(F); 

% STEP ONE: Perform wavelet transform on the signal being used for background power estimates
% NB: Background signal (bgsignal) is the time series (either row or column vector) used to estimate the background spectrum.

[B,T]=BOSC_tf(bgsignal,F,Fsample,width); % Compute the time-frequency (wavelet) spectrogram


% STEP TWO: Fit the background spectrum with a linear regression in log space to estimate mean power X frequency representation.
% Note that we use a shoulder based on the lowest frequency (longest period)

bgshoulder=ceil(width*Fsample/min(F)); % round up with ceil.m to be safe

bgfitMethod="optimized"; % choose fit method, options are: "standard", "robustfit", "median", "highpower", "freqsubset", "optimized"
[pv,meanpower]=BOSC_bgfit(F,B(:,(bgshoulder+1):(end-bgshoulder)),bgfitMethod);

% NB: meanpower is the estimated background power spectrum (geometric mean power as a function of frequency)
% NB: pv contains the slope and y-intercept of the regression line

% STEP THREE: Calculate the threshold values to use for detection

[powthresh,durthresh]=BOSC_thresholds(Fsample,percentilethresh,numcyclesthresh,F,meanpower);

% *** Hint: At this stage, it is a good idea to cross-check the background power spectrum fit (see PLOT #1: Power spectrum and background spectrum fit)

ks_d=BOSC_compare_chi2(F,B(:,(bgshoulder+1):(end-bgshoulder)),meanpower); % (optional) ks test to check background fit quality

% STEP FOUR: Set the target signal in which oscillations will be detected.
% The variable called "eegsignal" should contain more EEG signal than you want to analyze. starttime and endtime should be set to mark the bounds of a single trial of interest within eegsignal. Make sure there is enough additional signal before and after starttime and endtime (within eegsignal) to be able to include the shoulder to avoid edge artifacts (see README)
% NB: A safe way to calculate the shoulder:

shoulder=ceil((width+numcyclesthresh)*Fsample./F);

eegsignal=bgsignal; % eegsignal is assumed to include the time segment of interest, plus additional signal on either side.
starttime=max(shoulder)+1; % starttime contains the sample at which the trial of interest starts
endtime=length(bgsignal)-max(shoulder); % endtime contains the sample at which the trial of interest ends

DETECTED=[]; % reset matrix between uses
for f=1:length(F) % Loop through all frequencies (if all frequencies are desired; otherwise, specify a subset of frequencies)
  targetsignal=eegsignal((starttime-shoulder(f)):(endtime+shoulder(f)));
  % compute the time-frequency (wavelet) spectrogram, aka "scalogram"
  [Btarget,Ttarget]=BOSC_tf(targetsignal,F(f),Fsample,width); 

  % detect oscillations at frequency F(f)
  detected=BOSC_detect(Btarget,powthresh(f),durthresh(f),Fsample);
  detected=detected((shoulder(f)+1):(end-shoulder(f))); % *** Important: strip off the shoulders
  DETECTED(f,:)=detected; % accumulate all frequencies
end % frequency loop

Pepisode=mean(DETECTED,2); % Pepisode as a function of frequency. This is a useful summary measure.

% PLOT #1: Power spectrum and background spectrum fit

xf=1:4:length(F); % to get log-sampled frequency tick marks
plot(1:length(F),mean(log10(B(:,(bgshoulder+1):(end-bgshoulder))),2),'ko-',1:length(F),log10(meanpower),'r');
set(gca,'XTick',xf,'XTickLabel',F(1:4:end));
ylabel('Log(Power) [dB]');
xlabel('Frequency [Hz]');
xlim([0,22]);
%}

% PLOT #2: Pepisode(f) for a detected segment
%{
bar(1:length(F),Pepisode);
set(gca,'XTick',xf,'XTickLabel',F(1:4:end));
ylabel('P_{episode}'); xlabel('Frequency [Hz]');
%}

% PLOT #3: Raw-trace with superimposed BOSC-detected oscillations at a frequency of interest
%{
finterest=14; % Define a frequency of interest by indexing into F (i.e., not in Hz)
targetsignal=eegsignal(starttime:endtime);
% create copy of target data and find the detected values that = 1
osc=targetsignal; osc(find(DETECTED(finterest,:)~=1))=NaN; 			
t=(starttime:endtime)/Fsample;
h=plot(t,targetsignal,'b',t,osc,'g'); 
% Plot and label figure
ylabel('Voltage [\mu V]'); xlabel('Time [s]');
set(h(2),'LineWidth',2, 'Color', 'r'); % Colour detected oscillations (at frequency F(finterest)) red and use a thicker line
set(h(1),'LineWidth',1,'Color', 'k');  % Draw remaining signal in a thinner, black line
title(sprintf('P_{episode}(%.2f Hz)=%.2f',F(finterest),Pepisode(finterest))); % Optionally, put the frequency of interest and Pepisode in the title
axis tight;
%}

% PLOT #4: Time-frequency plot of BOSC-detected oscillations
%{
% This plots a black background and white denotes detected oscillations
% To reverse simply invert the values in the matrix: DETECTED=1-DETECTED;
imagesc(t,1:size(DETECTED,1),DETECTED);
colormap(gray);
set(gca,'YTick',xf,'YTickLabel',F(xf));
ylabel('Frequency [Hz]');
xlabel ('Time [s]');
axis xy;
%}

%{
% pdf and cdf of frequency
Bsorted=sort(B(:,(bgshoulder+1):(end-bgshoulder)),2);
f=16; % set frequency index
upperlim=round(chi2inv(0.99,2)); % approximate upper bound for x-axis based on chi-square percentile
% plot pdf of frequency
subplot(1,2,1);
hold on;
histogram(Bsorted(f,:),[0:1:upperlim,inf]*meanpower(f)/exp(psi(1)+log(2)),'Normalization','probability','FaceAlpha',1);
plot(Bsorted(f,:),pdf(makedist('Gamma','a',1,'b',2),Bsorted(f,:)/meanpower(f)*exp(psi(1)+log(2))),'r');
hold off;
xlim([0,(upperlim+1)*meanpower(f)/exp(psi(1)+log(2))]);
ylabel('Probalility');
xlabel('Power [\muV^2/Hz]');
% plot cdf of frequency
subplot(1,2,2);
p=plot(Bsorted(f,:),((1:length(Bsorted(f,:)))-0.5)/(length(Bsorted(f,:))),Bsorted(f,:),cdf(makedist('Gamma','a',1,'b',2),Bsorted(f,:)/meanpower(f)*exp(psi(1)+log(2))),'r');
xlim([0,upperlim*meanpower(f)/exp(psi(1)+log(2))]);
ylabel('Probalility');
xlabel('Power [\muV^2/Hz]');
legend([p(1),p(2)],'empirical','theoretical','Location','southeast');
%}