function y = dobp1Filter2(x, Fs, fcl, fch)
%DOFILTER Filters input x and returns output y.
%band pass filter designed with cut off frrequencies of 0.5 hz and 350 hz
%needs to be called by the filterbp function
% fcl=150;  %% lower cutoff frequency in Hz
% fch=500;  %% upper cutoff frequency in Hz

[B_low,A_low] = butter(6,2*[fcl fch]./Fs,'bandpass');
y = filtfilt(B_low,A_low,x);