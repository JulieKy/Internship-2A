function y = dobp1Filter(x, Fs)
%DOFILTER Filters input x and returns output y.
%band pass filter designed with cut off frrequencies of 0.5 hz and 350 hz
%needs to be called by the filterbp function

persistent Hd;

if isempty(Hd)
    
%     % The following code was used to design the filter coefficients:
%     %
%     N     = 4;      % Order
%     F3dB1 = 0.5;    % First
%     F3dB2 = 350;    % Second
%     % Fs    = 48000;  % Sampling Frequency
%     
%     h = fdesign.bandpass('n,f3db1,f3db2', N, F3dB1, F3dB2, Fs/2);
%     
%     Hd = design(h, 'butter', ...
%         'SystemObject', true);
    
    Hd = dsp.BiquadFilter( ...
        'Structure', 'Direct form II', ...
        'SOSMatrix', [1 0 -1 1 -1.99990744066141 0.999907444956381; 1 0 -1 1 ...
        -1.93540868194071 0.937435881661779], ...
        'ScaleValues', [0.0225116271536921; 0.0225116271536921; 1]);
end

s = double(x);
y = step(Hd,s);



