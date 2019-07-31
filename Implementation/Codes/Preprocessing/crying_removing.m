function [xsc] = crying_removing(xs)
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
[labels, coef_KAPPA]=labelling(observators,samples, end_sample, window, overlap);

%% CS POWER RATIO

%% NCS POWER RATIO

xsc=xs; % Need to be changed!
end

