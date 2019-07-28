%% PHASE 1

filename = uigetfile('*.m4a');
[y, Fs] = audioread(filename);
info = audioinfo(filename);
%sound(y,Fs)
if Fs > 16000
    yResamp = resample(y, 16000, Fs);
end

%get time of the audio signal
t =1/Fs:1/Fs:size(y)/Fs;
tResamp =1/16000:1/16000:size(yResamp)/16000;

% yCos = cos(2*pi*1000*tResamp);
% figure
% plot(tResamp(1:5000), yCos(1:5000))

% plot resampled audio file
figure("Name","Orignial Input","NumberTitle","off")
plot(tResamp, yResamp)
xlabel('Time(s)')
ylabel('Amplitude')

%% DEFINITIONS OF SOME VARIABLES

% define channel bands and other variables
channel = {{100 250},{250 500}, {500 1000}, {1000 2000}, {2000 4000}, {4000 8000}}; 
% aPassband = {{250, 500}};
lowestPassband = {{100, 250}};
% highestPassband = {{4000, 8000}};
% filterBank = {@WindowFIR, @gEquiFIR, @Chebyshev1IIR, @ButterWorthIIR};

%% FILTER

% collecting filtered audio, specify second argument for passband
% uncomment any section to filter and plot
butterWorthIIRChannels = filterWith(@ButterWorthIIR, channel, yResamp);
% plotFilteredAudio(butterWorthIIRChannels, tResamp)
% windowFIRChannels = filterWith(@WindowFIR, channel, yResamp);
% plotFilteredAudio(windowFIRChannels, tResamp)
% gEquiFIRChannels = filterWith(@gEquiFIR, lowestPassband, yResamp);
% plotFilteredAudio(gEquiFIRChannels, tResamp)
% chebyshevlIIRChannels = filterWith(@Chebyshev1IIR, lowestPassband, yResamp);
% plotFilteredAudio(chebyshevlIIRChannels, tResamp)

%% RECTIFY

RECbutterWorthIIRChannels = rectify(butterWorthIIRChannels);
% plotFilteredAudio(RECwindowFIRChannels, tResamp)
% RECwindowFIRChannels = rectify(windowFIRChannels);
% plotRectifiedAudio(RECwindowFIRChannels, tResamp)
% RECgEquiFIRChannels = rectify(windowFIRChannels);
% plotFilteredAudio(RECgEquiFIRChannels, tResamp)
% RECChebyshev1Channels = rectify(Chebyshev1Channels);
% plotFilteredAudio(RECChebyshev1Channels, tResamp)

%% LOWPASS/ENVELOPING

ENVbutterWorthIIRChannels = envelope(RECbutterWorthIIRChannels, WindowLowPass);
% plotEnvelopedSignal(ENVbutterWorthIIRChannels, tResamp);
% ENVwindowFIRChannels = envelope(RECwindowFIRChannels,WindowLowPass);
% plotEnvelopedSignal(ENVwindowFIRChannels, tResamp);
% ENVgEquiFIRChannels = envelope(RECgEquiFIRChannels, WindowLowPass);
% plotEnvelopedSignal(ENVgEquiFIRChannels, tResamp);
% ENVChebyshevlFIRChannels = envelope(RECChebyshev1Channels, WindowLowPass);
% plotEnvelopedSignal(ENVChebyshevlFIRChannels, tResamp);

%% Phase 3

centerFreq = getCenterFrequency(channel);
cosFunc = getCos(centerFreq, tResamp);
modulatedSig = ampModulation(cosFunc, ENVbutterWorthIIRChannels);
output = superposition(modulatedSig);
plotCos(output, tResamp)
audiowrite('finalButter.wav', output, 16000);
audiowrite('finalWindow.wav', output, 16000);
audiowrite('finalEqui.wav', output, 16000);
audiowrite('finalChebyshev.wav', output, 16000);

%% FUNCTIONS

function modulatedSig=ampModulation(cosSig, envSig)
modulatedSig={};
for channelIndex=1:size(envSig, 2)
    modulatedSig = [modulatedSig, times(envSig{channelIndex}, cosSig{channelIndex}')];
end
end

function allChannelSig=superposition(modulatedSigs)
allChannelSig = zeros(1, length(modulatedSigs{1}));
% for channelIndex=1:size(modulatedSigs,2)
%     allChannelSig = allChannelSig'+modulatedSigs{channelIndex};
%     disp(channelIndex)
allChannelSig = allChannelSig' + modulatedSigs{1};
allChannelSig = allChannelSig + modulatedSigs{2};
allChannelSig = allChannelSig + modulatedSigs{3};
allChannelSig = allChannelSig + modulatedSigs{4};
allChannelSig = allChannelSig + modulatedSigs{5};
allChannelSig = allChannelSig + modulatedSigs{6};
% end
end

function plotCos(signal, tResamp)
%     twoPeriods = 2*(1/centerFreq{1});
% for channelIndex=1:size(signal, 2)
%     figure
%     plot(tResamp, signal{channelIndex})
% end
    figure
    plot(tResamp, signal)
end

function cosFunc=getCos(centerFreq, tResamp)
cosFunc = {};
for channel=1:length(centerFreq)
    cosFunc = [cosFunc, cos(2*pi*centerFreq(channel)*tResamp)];
end
end

function centerFreq=getCenterFrequency(cutOffFreqs)
centerFreq = [];
for range =1:length(cutOffFreqs)
    passband = cutOffFreqs{range};
    centerFreq = [centerFreq, sqrt(passband{1}*passband{2})];
end
end

% Enveloping signal
function envelopedAudio=envelope(RECchannelResults, lowpassFilter)
envelopedAudio = {};
for channelIndex=1:length(RECchannelResults)
   envelopedAudio =[envelopedAudio, {filter(lowpassFilter, RECchannelResults{channelIndex})}];
end
end

% Plotting enveloped signal
function plotEnvelopedSignal(ENVchannelResults, tResamp)
for channelIndex=1:length(ENVchannelResults)
    figure('Name', 'Enveloped Signal')
    plot(tResamp,ENVchannelResults{channelIndex})
    xlabel('Time(s)')
    ylabel('Amplitude')
end  
end

% Rectifying signal
function rectifiedAudio=rectify(channelResults)
rectifiedAudio = {};
for channelIndex=1:size(channelResults,1)
    if size(channelResults,1) ~= 1
        channelDetails = channelResults(channelIndex, 1:end);
    else
        channelDetails = channelResults(1:end);
    end
    rectifiedAudio = [rectifiedAudio, abs(channelDetails{2})];
end
end

% plotting rectified signal
function plotRectifiedAudio(RECchannelResults, tResamp)
for channelIndex=1:length(RECchannelResults)
    figure('Name','Rectified Signal')
    plot(tResamp,RECchannelResults{channelIndex})
    xlabel('Time(s)')
    ylabel('Amplitude')
end  
end

% filter the audio signal with specified filter and passband range
function out = filterWith(filterName,channel,y)
if nargin<3
    channel = {{100 250},{250 500}, {500 1000}, {1000 2000}, {2000 4000}, {4000 8000}};
end
out = {};
for range =1:length(channel)
    tic
    filterOutput = filter(filterName(channel{range}), y);
    toc
%     out = [out; {channel{range},filter(filterName(channel{range}),y)}];
    out = [out; {channel{range},filterOutput}];
end
end

% plotting the filtered channel signals
function plotFilteredAudio(channelResults, tResamp)
for channelIndex=1:size(channelResults,1)
    if size(channelResults,1) ~= 1
        channelDetails = channelResults(channelIndex, 1:end);
    else
        channelDetails = channelResults(1:end);
    end
    channelPassband = channelDetails{1,1};
    figure("NumberTitle","off","Name",strcat("Channel",string(channelPassband).join))
    plot(tResamp, channelDetails{1,2});
    xlabel('Time(s)')
    ylabel('Amplitude')
end
end

% Lowpass filter for enveloping
function Hd = WindowLowPass
%WINDOWLOWPASS1 Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.4 and Signal Processing Toolbox 8.0.
% Generated on: 15-Jul-2019 17:21:56

% FIR Window Lowpass filter designed using the FIR1 function.

% All frequency values are in Hz.
Fs = 16000;  % Sampling Frequency

N    = 50;      % Order
Fc   = 400;      % Cutoff Frequency
flag = 'scale';  % Sampling Flag
Beta = 0.5;      % Window Parameter

% Create the window vector for the design algorithm.
win = kaiser(N+1, Beta);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, Fc/(Fs/2), 'low', win, flag);
Hd = dfilt.dffir(b);

% [EOF]
end

% Bank of Filters are defined from here to end of file
function Hd = WindowFIR(Fc)
%WINDOW Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.4 and DSP System Toolbox 9.6.
% Generated on: 14-Jul-2019 15:49:28

% FIR Window Bandpass filter designed using the FIR1 function.

% All frequency values are in Hz.
Fs = 16000;  % Sampling Frequency

N    = 50;      % Order
Fc1  = Fc{1};      % First Cutoff Frequency
Fc2  = Fc{2};      % Second Cutoff Frequency
flag = 'scale';  % Sampling Flag
Beta = 0.5;      % Window Parameter
% Create the window vector for the design algorithm.
win = kaiser(N+1, Beta);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, [Fc1 Fc2-1]/(Fs/2), 'bandpass', win, flag);
Hd = dfilt.dffir(b);

% [EOF]
end

function Hd = gEquiFIR(Fc)
%GEQUI Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.4 and DSP System Toolbox 9.6.
% Generated on: 14-Jul-2019 15:25:47

% Generalized REMEZ FIR Bandpass filter designed using the FIRGR function.

% All frequency values are in Hz.
Fs = 16000;  % Sampling Frequency

Fstop1 = Fc{1}-50;             % First Stopband Frequency
Fpass1 = Fc{1};             % First Passband Frequency
Fpass2 = Fc{2};             % Second Passband Frequency
Fstop2 = Fc{2}+50;             % Second Stopband Frequency
Dstop1 = 0.0001;          % First Stopband Attenuation
Dpass  = 0.057501127785;  % Passband Ripple
Dstop2 = 0.0001;          % Second Stopband Attenuation
dens   = 200;             % Density Factor

% Calculate the coefficients using the FIRGR function.
b  = firgr('mineven', [0 Fstop1 Fpass1 Fpass2 Fstop2 Fs/2]/(Fs/2), [0 ...
           0 1 1 0 0], [Dstop1 Dpass Dstop2], {dens});
Hd = dfilt.dffir(b);

% [EOF]

end

function Hd = Chebyshev1IIR(Fc)
%CHEBYSHEV1 Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.4 and DSP System Toolbox 9.6.
% Generated on: 14-Jul-2019 15:53:27

% Chebyshev Type I Bandpass filter designed using FDESIGN.BANDPASS.

% All frequency values are in Hz.
Fs = 1600;  % Sampling Frequency

N      = 50;   % Order
Fpass1 = Fc{1};  % First Passband Frequency
Fpass2 = Fc{2};  % Second Passband Frequency
Apass  = 1;    % Passband Ripple (dB)

% Construct an FDESIGN object and call its CHEBY1 method.
h  = fdesign.bandpass('N,Fp1,Fp2,Ap', N, Fpass1, Fpass2, Apass, Fs);
Hd = design(h, 'cheby1');

% [EOF]
end

function Hd = ButterWorthIIR(Fc)
%Fc should be a list of size 2
%TASK4IIR Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.4 and Signal Processing Toolbox 8.0.
% Generated on: 12-Jul-2019 16:27:16

% Butterworth Bandpass filter designed using FDESIGN.BANDPASS.

% All frequency values are in Hz.
Fs = 16000;  % Sampling Frequency
Fc1 = Fc{1};    % First passband Frequency
Fc2 = Fc{2};    % Second passband Frequency
N   = 50;   % Order

% Construct an FDESIGN object and call its BUTTER method.
% Fc(1) defines low cut off Fc(2) defines high cut off
h  = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, Fs);
Hd = design(h, 'butter');

% [EOF]

end
