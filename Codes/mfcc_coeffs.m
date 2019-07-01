function [output_mean_mfcc] = mfcc_coeffs(y, fn, sample)
%%Gives the mean mfcc coefficients of the first 6 mfcc coefficients


%% INPUT AND OUTPUT
% any input audio signal
% can be filtered or raw

% input = x 'audio signal'
% Fs 'sampling frequency
% output =  mean mfcc coefficients of the first 6 mfcc coefficients


%% INITIALISATION

total_samples = length(y);
numberCoeffs = 6;


%% COMPUTATION OF MFCC COEFFICIENTS

% finding 6 mfcc coefficients for full signal/the input signal
coeffs = mfcc(y,fn, 'NumCoeffs', numberCoeffs);

%% RESULT

% output of the average first 6 MFCC (coeff 1-2-3-4-5-6) 
 output_mean_mfcc = mean(coeffs(:,1:numberCoeffs));
