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
window=3;
overlap=25/100;

% -- Parameters for the power ratio
pass_band=[0:2000];
band_width=100;


path = pwd;

%% LABELLING
[label_final, coef_KAPPA]=labelling(observators,samples, end_sample, window, overlap);


%% POWER RATIO

pxx_NCS=[];
band_NCS=[];
PR_NCS=[];

pxx_CS=[];
band_CS=[];
PR_CS=[];

% For every signal
for i = 1:length(names_cell)
    
    % -- Reading the signals
    tempName=names_cell{i};
    [x,Fs]= audioread([path,'\..\Data\Samples_Belle\',tempName]);
    
    % -- Resampling to 4000 Hz
    xs=resample(x,4000,Fs);
    fn=4000;
    
    % -- Shorten to 60s
    time_sample=60;
    xss=xs(1:time_sample*fn,1);
    
    % -- NCS power ratio
    flag_section=0; % 0 for NCS
    [pxx_NCS_signal, band_NCS_signal, PR_NCS_signal, freq]=power_ratio_band(xss, i, fn, window, overlap, label_final, pass_band, band_width, flag_section);
    pxx_NCS=[pxx_NCS_signal; pxx_NCS];
    band_NCS=[band_NCS_signal; band_NCS];
    PR_NCS=[PR_NCS_signal; PR_NCS];
    f=freq;
    
    % -- CS power ratio
    flag_section=1; % 1 for CS
    [pxx_CS_signal, band_CS_signal, PR_CS_signal, freq]=power_ratio_band(xss, i, fn, window, overlap, label_final, pass_band, band_width, flag_section);
    if (pxx_CS_signal~=0) % CS present in the signal
        pxx_CS=[pxx_CS_signal; pxx_CS];
        band_CS=[band_CS_signal; band_CS];
        PR_CS=[PR_CS_signal; PR_CS];
    end
    
end

% NCS: average on all signals
pxx_NCS_mean=mean(pxx_NCS);
band_NCS_mean=mean(band_NCS);
PR_NCS_mean=mean(PR_NCS);

% CS: average on all signals
pxx_CS_mean=mean(pxx_CS);
band_CS_mean=mean(band_CS);
PR_CS_mean=mean(PR_CS);


%% CROSS VALIDATION


%% DISPLAY

% Display annoted labels NCS and CS
signal_n=7;
display_NCS_CS_annotations(signal_n,label_final, window, overlap)

% %  Display the periodograms of annotated NCS and CS
% display_PR_NCS_CS_interquartiles(f,pxx_NCS, pxx_CS, pxx_NCS_mean, pxx_CS_mean, band_width, pass_band, band_NCS_mean, band_CS_mean);


threshold=1;
end
