function simsignal = simulate(F)
%% simulation parameters

cfg.simParams.rhythmFreq    = [10 20 16];         % oscillatory frequencies
cfg.simParams.amplitude     = [0 2 4 6 8 12 16];     % simulated signal SNR
cfg.simParams.cycles        = [600 1200 960];     % simulated signal durations [in cycles]
cfg.simParams.segmentDur    = 270;                % duration of total segement [in seconds]
cfg.simParams.fsample       = 500;                % sampling rate of simulated signal [in Hz]

%% create background data
%   
% For reproducibility, we use a fixed seed.
randn('seed',13);
% generate 1/f background with the function 'f_alpha_gaussian' from the CNOISE toolbox
% https://people.sc.fsu.edu/~jburkardt/m_src/cnoise/cnoise.html
bckgrnd_filt = f_alpha_gaussian(cfg.simParams.segmentDur*cfg.simParams.fsample,1,1)';
% bandpass filter signal (consecutive low + high-pass filter)
[B,A]  = butter(4,70/(cfg.simParams.fsample/2),'low'); 
bckgrnd_filt = filtfilt(B,A,bckgrnd_filt); clear A B
[B,A]  = butter(4,.5/(cfg.simParams.fsample/2),'high'); 
bckgrnd_filt = filtfilt(B,A,bckgrnd_filt); clear A B
%% create simulated data structure with specified power and duration
simsignal=[];
%for a = 1:length(cfg.simParams.amplitude)
for a=4 % potential for loop through amplitude
    data.time = [(1/cfg.simParams.fsample):(1/cfg.simParams.fsample):cfg.simParams.segmentDur];
    % generate alpha in the middle of the segment
    rhythmTime = round((cfg.simParams.cycles(1)./cfg.simParams.rhythmFreq(1)),3);
    % simulate rhythms as symmetrical around the center
    timeNew = round(rhythmTime./(1/cfg.simParams.fsample),0);
    if mod(timeNew,2) ~= 0
    timeNew = timeNew + 1;
    rhythmTime = timeNew.*(1/cfg.simParams.fsample);
    else rhythmTime = timeNew.*(1/cfg.simParams.fsample);
    end; clear timeNew;
    rhythmTimeVector = [(1/cfg.simParams.fsample):(1/cfg.simParams.fsample):rhythmTime];
    % oscillatory period from 80-140 s, centered at 110 s
    rhythmIdxVector1 = (110*cfg.simParams.fsample)-(numel(rhythmTimeVector)/2)+1:(110*cfg.simParams.fsample)+(numel(rhythmTimeVector)/2);
    % oscillatory period from 200-260 s, centered at 230 s
    rhythmIdxVector2 = (230*cfg.simParams.fsample)-(numel(rhythmTimeVector)/2)+1:(230*cfg.simParams.fsample)+(numel(rhythmTimeVector)/2);
    rhythmIdxVector=[rhythmIdxVector1,rhythmIdxVector2]; % combining both together
    % filter entire signal between 8 and 12 Hz (6th order butterworth) (not locally on alpha)
    [tmp_b,tmp_a] = butter(6, [8, 12]/(250/2), 'bandpass');
    tmp_bpsignal = filter(tmp_b,tmp_a,squeeze(bckgrnd_filt)); clear tmp_a tmp_b;
    % scale rhythm by noise background power
    amplitudeFromRMS = (sqrt(cfg.simParams.amplitude(a)*var(tmp_bpsignal))*sqrt(2)); clear tmp_bpsignal;
    % single peak
    simulatedRhythm = [sin(rhythmTimeVector*2*pi*cfg.simParams.rhythmFreq(1))*amplitudeFromRMS,sin(rhythmTimeVector*2*pi*cfg.simParams.rhythmFreq(1))*amplitudeFromRMS];
    % possible second peak
    %simulatedRhythm=simulatedRhythm+[sin(rhythmTimeVector*2*pi*cfg.simParams.rhythmFreq(2))*amplitudeFromRMS,sin(rhythmTimeVector*2*pi*cfg.simParams.rhythmFreq(2))*amplitudeFromRMS];
    % possible third peak, half the size
    %simulatedRhythm=simulatedRhythm+[sin(rhythmTimeVector*2*pi*cfg.simParams.rhythmFreq(3))*amplitudeFromRMS/2,sin(rhythmTimeVector*2*pi*cfg.simParams.rhythmFreq(3))*amplitudeFromRMS/2];

    simulatedRhythm_complete = zeros(1,numel(data.time));
    simulatedRhythm_complete(1,rhythmIdxVector) = simulatedRhythm;
    %{
    % simulating high power vals at low freqs (1-4 Hz)
    simhighPower=zeros(1,numel(data.time));
    rng(20,"twister"); % consistent rng seed
    highpwrRMSamplitude=2; % high amplitude (~cfg.simParams.amplitude=16) 
    for h=1:9 % 9 randomly placed high power rhythms, 15 s each
    startidx=randi(numel(data.time)-(15*cfg.simParams.fsample)); 
    simhighpowerRhythms=sin([startidx:(1/cfg.simParams.fsample):startidx+15]*2*pi*F(h))*highpwrRMSamplitude;
    simhighPower(startidx:(startidx+15*cfg.simParams.fsample))=simhighpowerRhythms; % put together as one signal
    end
    simulatedRhythm_complete=simulatedRhythm_complete+simhighPower; % add to other rhythms
    %}
    % effective segment = BG + rhythm %
    simsignal = [simsignal; bckgrnd_filt + simulatedRhythm_complete];
    clear simulatedRhythm_complete simulatedRhythm amplitudeFromRMS
end 
clear count rhythmTimeVector rhythmIdxVector rhythmTime a c k bckgrnd_filt

