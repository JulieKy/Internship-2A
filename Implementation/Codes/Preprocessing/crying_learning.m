function [threshold] = crying_learning(names_cell)
%CRYING_LEARNING:  Label the crying section, and learn where are the CS by using the Power Ratio Tool

%% INPUTS AND OUTPUTS
%  -- Inputs --
% names_cell:
% -- Outputs --


%% INITIALISATION
% -- Data
observators=3;
samples=37;
end_sample=60; % End of the signal (hypotesis: length of the signal=60s)

% -- Parameters for labelling
window=1;
overlap=0;

% -- Parameters for the power ratio
pass_band=[0:2000];
band_width=100;


path = pwd;

%% LABELLING
[label_final, coef_KAPPA]=labelling(observators,samples, end_sample, window, overlap);


%% POWER RATIO

% Initialisation
pxx_NCS=[]; band_NCS=[]; PR_NCS=[]; p25_NCS=[]; p75_NCS=[];
pxx_CS=[]; band_CS=[]; PR_CS=[]; p25_CS=[]; p75_CS=[];

% For every signal
for i = 1:length(names_cell)
    
    % -- Reading the signals
    tempName=names_cell{i};
    [x,Fs]= audioread([path,'\..\Data\Samples_Belle\',tempName]);
    disp('READ');
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
    
    % -- NCS power ratio
    flag_section=0; % 0 for NCS
    [pxx_NCS_signal, band_NCS_signal, PR_NCS_signal, freq, p25_NCS_signal, p75_NCS_signal]=power_ratio_band(xss, signal_n, fn, window, overlap, label_final, pass_band, band_width, flag_section);
    pxx_NCS=[pxx_NCS_signal; pxx_NCS];
    band_NCS=[band_NCS_signal; band_NCS];
    PR_NCS=[PR_NCS_signal; PR_NCS];
    f=freq;
    p25_NCS=[p25_NCS, p25_NCS_signal];
    p75_NCS=[p75_NCS, p75_NCS_signal];
    
    % -- CS power ratio
    flag_section=1; % 1 for CS
    [pxx_CS_signal, band_CS_signal, PR_CS_signal, freq, p25_CS_signal, p75_CS_signal]=power_ratio_band(xss, signal_n, fn, window, overlap, label_final, pass_band, band_width, flag_section);
    if (pxx_CS_signal~=0) % CS present in the signal
        pxx_CS=[pxx_CS_signal; pxx_CS];
        band_CS=[band_CS_signal; band_CS];
        PR_CS=[PR_CS_signal; PR_CS];
        p25_CS=[p25_CS, p25_CS_signal];
        p75_CS=[p75_CS, p75_CS_signal];
    end
    
end

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


%% DISPLAY

% Display annoted labels NCS and CS
signal_n=22;
% display_NCS_CS_annotations(signal_n,label_final, window, overlap)

%  Display the periodograms of annotated NCS and CS
display_PR_NCS_CS_interquartiles(f,pxx_NCS, pxx_CS, pxx_NCS_mean, pxx_CS_mean, band_width, pass_band, band_NCS_mean, band_CS_mean);


threshold=1;

end
