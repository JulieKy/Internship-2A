function [xsc] = crying_removing(xss, fn)
%CRYING_REMOVING:  Label the crying section, detect them and remove them

%% INPUTS AND OUTPUTS
%  -- Inputs --
% xs: input signal
% -- Outputs --
% xsc: input signal without crying sections

%% INITIALISATION
% -- Data
observators=2;
samples=37;
end_sample=60; % End of the signal (hypotesis: length of the signal=60s)

% -- Parameters for labelling
window=3;
overlap=25/100;

%% LABELLING
[label_final, coef_KAPPA]=labelling(observators,samples, end_sample, window, overlap);

%% POWER RATIO
pass_band=[0:1000];
band_width=100;

% -- For NCS
flag_section=0; % 0 for NCS
[pxx_mean_NCS, band_mean_NCS, PR_mean_NCS]=power_ratio_band(xss, fn, window, overlap, label_final, pass_band, band_width, flag_section);

% -- For CS
flag_section=1; % 1 for CS
[pxx_mean_CS, band_mean_CS, PR_mean_CS]=power_ratio_band(xss, fn, window, overlap, label_final, pass_band, band_width, flag_section);

xsc=xss; % Need to be changed!
end

