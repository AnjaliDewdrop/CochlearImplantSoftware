function Hd = WindowFIR
%WINDOW Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.4 and DSP System Toolbox 9.6.
% Generated on: 14-Jul-2019 15:49:28

% FIR Window Bandpass filter designed using the FIR1 function.

% All frequency values are in Hz.
Fs = 16000;  % Sampling Frequency

N    = 200;      % Order
Fc1  = 250;      % First Cutoff Frequency
Fc2  = 500;      % Second Cutoff Frequency
flag = 'scale';  % Sampling Flag
Beta = 0.5;      % Window Parameter
% Create the window vector for the design algorithm.
win = kaiser(N+1, Beta);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, [Fc1 Fc2]/(Fs/2), 'bandpass', win, flag);
Hd = dfilt.dffir(b);

% [EOF]
end
