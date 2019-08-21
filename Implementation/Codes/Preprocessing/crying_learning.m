function [final_threshold, band, label_final, window, overlap] = crying_learning(names_cell)
%CRYING_LEARNING:  Label the crying section, and learn where are the CS by using the Power Ratio Tool

%% INPUTS AND OUTPUTS
%  -- Inputs --
% names_cell: Name of the samples 
% -- Outputs --
% final_threshold: Threshold used for distinguish CS and NCS
% band: Frequency band used to compare the power ratio of CS and NCS
% label_final: Annotated labels of the signals
% window: Window used for labelling (annotations)
% overlap: Overlap used for labelling (annotations)

%% INITIALISATION
% -- Data
observators=3;
samples=37;
end_sample=60; % End of the signal (hypotesis: length of the signal=60s)

% -- Parameters for labelling
window_labelling=1;
overlap_labelling=0;

% -- Parameters for training
window_training=3;

% -- Parameters for the power ratio
pass_band=[0:2000];
band_width=100;

% -- Path
path = pwd;

% -- Initialisation training segments 
CS_pure=[]; NCS_pure=[]; CS_epsi=[]; NCS_epsi=[];

%% LABELLING
[label_final, coef_KAPPA]=labelling(observators,samples, end_sample, window_labelling, overlap_labelling);


%% TRAINING

% Initialisation
pxx_NCS=[]; band_NCS=[]; PR_NCS=[]; p25_NCS=[]; p75_NCS=[];
pxx_CS=[]; band_CS=[]; PR_CS=[]; p25_CS=[]; p75_CS=[];

%% -- Reading files and first preprocessing
% For every signal
for i = 1:length(names_cell)
    
    % -- Reading the signals
    tempName=names_cell{i};
    [x,Fs]= audioread([path,'\..\Data\Samples_Belle\',tempName]);
    disp('READ - crying_learning.m');
    disp(tempName);
    
    % Get the number of the recording by removing the '.mp3'
    strMP3 = sprintf('%s',tempName);
    ind=strfind(strMP3,'.');
    signal_n = str2num(strMP3(1:ind-1));
    
    % -- Resampling to 4000 Hz
    xs=resample(x,4000,Fs);
    fn=4000;
    
    % -- Shorten to 60s
    time_sample=60;
    xss=xs(1:time_sample*fn,1);
    
    %% -- Finding pure CS and NCS
    label_final_xss=label_final(signal_n, :);
    epsilon=0;
    [CS_pure_xss, NCS_pure_xss, CS_epsi_xss, NCS_epsi_xss] = CS_NCS_pure_window(xss, fn, signal_n, label_final_xss, window_training, epsilon, window_labelling);
    CS_pure=[CS_pure_xss; CS_pure];
    NCS_pure=[NCS_pure_xss; NCS_pure];
    CS_epsi=[CS_epsi_xss; CS_epsi];
    NCS_epsi=[NCS_epsi_xss; NCS_epsi];
end 

    %% -- Power ratio
    % CS power ratio
    [pxx_CS, band_CS, PR_CS, freq, p25_CS, p75_CS]=power_ratio_band2(pass_band, fn, CS_pure);
    
    % NCS power ratio
    [pxx_NCS, band_NCS, PR_NCS, freq, p25_NCS, p75_NCS]=power_ratio_band2(pass_band, fn, NCS_pure);
    

% NCS: average on all signals
pxx_NCS_mean=mean(pxx_NCS(pxx_NCS~=0));
band_NCS_mean=mean(band_NCS(band_NCS~=0));
PR_NCS_mean=mean(PR_NCS(PR_NCS~=0));
p25_NCS_mean=mean(p25_NCS(p25_NCS~=0));
p75_NCS_mean=mean(p75_NCS(p75_NCS~=0));

% CS: average on all signals
pxx_CS_mean=mean(pxx_CS(pxx_CS~=0));
band_CS_mean=mean(band_CS(band_CS~=0));
PR_CS_mean=mean(PR_CS(PR_CS~=0));
p25_CS_mean=mean(p25_CS(p25_CS~=0));
p75_CS_mean=mean(p75_CS(p75_CS~=0));

%% SPECTROGRAM
signal_n=15; % Signal wanted
wind_time_spec=0.5; % Window of 1 second
overlap_spec=0.25; % 25% overlap
start_time=0; end_time=15; % Part of the signal wanted
[s, f, t] = spectrogram_CS(signal_n, wind_time_spec, overlap_spec, start_time, end_time);


%% NCS&CS FEATURES EXTRACTION
% [zrc_CS, output_spectral_features_CS, periodogram_pks_features_CS, output_mean_mfcc_CS, output_lpc_CS, output_lsf_CS, zrc_NCS, output_spectral_features_NCS, periodogram_pks_features_NCS, output_mean_mfcc_NCS, output_lpc_NCS, output_lsf_NCS,] = NCS_CS_features_boxplot();

%% THRESHOLD DETERMINATION

% -- Powerband
band=[p25_CS_mean, p75_CS_mean];
powerband=[];

% For each signal
for i = 1:length(names_cell)
    
    % Reading the signals
    tempName=names_cell{i};
    [x,Fs]= audioread([path,'\..\Data\Samples_Belle\',tempName]);
    disp('READ - crying_learning.m THRESHOLD');
    disp(tempName);
    
    % Get the number of the recording by removing the '.mp3'
    strMP3 = sprintf('%s',tempName);
    ind=strfind(strMP3,'.');
    signal_n = str2num(strMP3(1:ind-1));
    
    % -- Resampling to 4000 Hz
    xs=resample(x,4000,Fs);
    fn=4000;
    
    % -- Shorten to 60s
    time_sample=60;
    xss=xs(1:time_sample*fn,1);
    
    % Find the band power of each CS and NCS of the signal
    labels_x= label_final(signal_n, :);
    powerband_signal = powerband_segment(band, fn, xss, labels_x, window, overlap);
    powerband=[powerband; powerband_signal];
end

% -- ROC
nb_thresholds=500; % Number of thresholds for the ROC
label_annotated=powerband(:,1);
powerband=powerband(:,2);
[fpr, tpr, final_threshold] = threshold_ROC( nb_thresholds, label_annotated, powerband); % Compute the ROC



%% DISPLAY

% Display annoted labels NCS and CS
signal_n=22;
display_NCS_CS_annotations(signal_n,label_final, window, overlap)

%  Display the periodograms of annotated NCS and CS
display_PR_NCS_CS_interquartiles(f,pxx_NCS, pxx_CS, pxx_NCS_mean, pxx_CS_mean, band_width, pass_band, band_NCS_mean, band_CS_mean);
end
