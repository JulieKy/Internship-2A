function [y] = filterbp(x,Fs)
%% Gives the filtered signal
% input = raw audio signal read as x
% Fs = sampling frequency
% output = band pass filtered signal with cut off frequencies of 0.5hz, and
% 350 hz

y = dobp1Filter2(x, Fs,fcl,fch);%calling band pass filtering function
audiowrite('filtered_signal.wav',y,Fs); %writing the filtered signal to an audio file

end

