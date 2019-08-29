function [output_mean_mfcc] = mfcc_coeffs(y, fn, label_training, length_labels_training, signal_n)
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
coeffs = mfcc(y',fn, 'NumCoeffs', numberCoeffs);

%% REMOVE THE COEEFICIENT CORRESPONDING TO CSs
labels_withPadding=label_training(signal_n, :);
length_label=length_labels_training(signal_n); 
labels=labels_withPadding(1:length_label);
NCS=1-labels;
for i=1:numberCoeffs
    coeff=coeffs(:,i)
    coeffs_NCS(:, i)=coeff((coeff.*NCS)~=0);
end 


%% RESULT

% output of the average first 6 MFCC (coeff 1-2-3-4-5-6) 
 output_mean_mfcc = mean(coeffs_NCS(:,1:numberCoeffs));