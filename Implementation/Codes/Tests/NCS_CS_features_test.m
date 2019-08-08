% Boxplot NCS/CS features

%% INITIALISATION PATH
clear all, close all, clc, dbstop if error;
addpath(genpath('..\..\')); % to have access to sample folder
init = 0; % optional, to generate a new Excel File
excelFile = 'NewFeaturesAnalysis_f'; % name of the Excel file to store features
path = pwd; % current path
data_dir=[path,'\..\..\Data\Samples_Belle\'];


%% READY TO READ FILES
dinfo = dir(data_dir);
names_cell1 = {dinfo.name};
% choose a valid file name
j=0;
for i=1:size(names_cell1,2)
    if length(names_cell1{i})>2
        j=j+1;
        names_cell{j}=names_cell1{i};
    end
end
lengthTot=j;

%% LABELLING NCS/CS
% -- Data
observators=3;
samples=37;
end_sample=60; % End of the signal (hypotesis: length of the signal=60s)

% -- Parameters for labelling
window=3;
overlap=25/100;

% -- Labelling
[label_final, coef_KAPPA]=labelling(observators,samples, end_sample, window, overlap);

% -- Initialisation for feature storage
zrc_CS=[]; output_spectral_features_CS=[]; periodogram_pks_features_CS=[]; output_mean_mfcc_CS=[]; output_lpc_CS=[]; output_lsf_CS=[];
zrc_NCS=[]; output_spectral_features_NCS=[]; periodogram_pks_features_NCS=[]; output_mean_mfcc_NCS=[]; output_lpc_NCS=[]; output_lsf_NCS=[];

%% FEATURES EXTRACTION

% -- For each recording
for i = 1:lengthTot
    
    %% -- Reading recordings
    close all;
    tempName=names_cell{i};
    
    % Get the number of the recording by removing the '.mp3'
    strMP3 = sprintf('%s',tempName);
    ind=strfind(strMP3,'.');
    signal_n = str2num(strMP3(1:ind-1));
    str=sprintf(' -- tempname: %s & signal_n: %d --\n', tempName, signal_n);
    disp(str)
    
    [x,Fs]= audioread([path,'\..\..\Data\Samples_Belle\',tempName]); % read current file
    
    %% -- Resampling to 4000 Hz
    xs=resample(x,4000,Fs);
    fn=4000;
    
    %% -- Shorten to 60s
    time_sample=60;
    xss=xs(1:time_sample*fn,1);
    
    % -- Variables
    N = length(xss);
    time_axis = (1:N)/fn;
    
    
    %% -- Finding CS on the signal
    flag_section=1; % CS
    locs=find(label_final(i,:)==flag_section); % Locations of CS
    
    if isempty(locs)==0 % There are CS on the signal
        
        % -- Start time of the labels (for each window)
        start_time=locs*(window-window*overlap)-1;
        
        % -- Start sample of the labels (for each window)
        start_sample=start_time*fn+1;
        label_duration=window*fn; % Number of samples in a window
        
        for n_section=1:length(locs)
            [xss_section,time_axis_section] = label2signal(xss, n_section, start_sample, label_duration, time_axis);
            
            %% -- Features extraction
            zrc = temporal_features(xss_section,fn, tempName); % Temporal features
            [output_spectral_features, periodogram_pks_features, pxx, f, foct, spower, I, S] = spectral_features(xss_section,fn); % See Fae's comment
            output_mean_mfcc = mfcc_coeffs(xss_section, fn); % MFCCs coefficient
            [output_lpc, output_lsf] = lpc_lsf_coeff(xss_section, fn); % LPC and LFC coefficient
            
            %% -- Storage of features
            % A new section = a new line
            zrc_CS=[zrc; zrc_CS];
            output_spectral_features_CS=[output_spectral_features'; output_spectral_features_CS];
            periodogram_pks_features_CS=[periodogram_pks_features'; periodogram_pks_features_CS];
            output_mean_mfcc_CS=[output_mean_mfcc; output_mean_mfcc_CS];
            output_lpc_CS=[output_lpc; output_lpc_CS];
            output_lsf_CS=[output_lsf; output_lsf_CS];
        end
    end
    meanPSD_CS=output_spectral_features_CS(:,1);
    stdPSD_CS=output_spectral_features_CS(:,2);
    medPSD_CS=output_spectral_features_CS(:,3);
    bw_CS=output_spectral_features_CS(:,4);
    p25_CS=output_spectral_features_CS(:,5);
    p75_CS=output_spectral_features_CS(:,6);
    IQR_CS=output_spectral_features_CS(:,7);
    TP_CS=output_spectral_features_CS(:,8);
    p100_200_CS=output_spectral_features_CS(:,9);
    p200_400_CS=output_spectral_features_CS(:,10);
    p400_800_CS=output_spectral_features_CS(:,11);
    spectrum_slope2_CS=output_spectral_features_CS(:,12);
    r_square2_CS=output_spectral_features_CS(:,13);
    
    
    %% -- Finding NCS on the signal
    flag_section=0; % NCS
    locs=find(label_final(i,:)==flag_section); % Locations of NCS
    
    if isempty(locs)==0 % There are NCS on the signal
        
        % -- Start time of the labels (for each window)
        start_time=locs*(window-window*overlap)-1;
        
        % -- Start sample of the labels (for each window)
        start_sample=start_time*fn+1;
        label_duration=window*fn; % Number of samples in a window
        
        for n_section=1:length(locs)
            [xss_section,time_axis_section] = label2signal(xss, n_section, start_sample, label_duration, time_axis);
            
            %% -- Features extraction
            zrc = temporal_features(xss_section,fn, tempName); % Temporal features
            [output_spectral_features, periodogram_pks_features, pxx, f, foct, spower, I, S] = spectral_features(xss_section,fn); % See Fae's comment
            output_mean_mfcc = mfcc_coeffs(xss_section, fn); % MFCCs coefficient
            [output_lpc, output_lsf] = lpc_lsf_coeff(xss_section, fn); % LPC and LFC coefficient
            
            %% -- Storage of features
            % A new section = a new line
            zrc_NCS=[zrc; zrc_NCS];
            output_spectral_features_NCS=[output_spectral_features'; output_spectral_features_NCS];
            periodogram_pks_features_NCS=[periodogram_pks_features'; periodogram_pks_features_NCS];
            output_mean_mfcc_NCS=[output_mean_mfcc; output_mean_mfcc_NCS];
            output_lpc_NCS=[output_lpc; output_lpc_NCS];
            output_lsf_NCS=[output_lsf; output_lsf_NCS];
        end
    end
end

    meanPSD_NCS=output_spectral_features_NCS(:,1);
    stdPSD_NCS=output_spectral_features_NCS(:,2);
    medPSD_NCS=output_spectral_features_NCS(:,3);
    bw_NCS=output_spectral_features_NCS(:,4);
    p25_NCS=output_spectral_features_NCS(:,5);
    p75_NCS=output_spectral_features_NCS(:,6);
    IQR_NCS=output_spectral_features_NCS(:,7);
    TP_NCS=output_spectral_features_NCS(:,8);
    p100_200_NCS=output_spectral_features_NCS(:,9);
    p200_400_NCS=output_spectral_features_NCS(:,10);
    p400_800_NCS=output_spectral_features_NCS(:,11);
    spectrum_slope2_NCS=output_spectral_features_NCS(:,12);
    r_square2_NCS=output_spectral_features_NCS(:,13);
    
    %% DISPLAY
    w=2;
    h=7;
    
    meanPSD=[meanPSD_CS;meanPSD_NCS]'
    stdPSD=[stdPSD_CS;stdPSD_NCS];
    medPSD=[medPSD_CS;medPSD_NCS];
    bw=[bw_CS;bw_NCS];
    p25=[p25_CS;p25_NCS];
    p75=[p75_CS;p75_NCS];
    IQR=[IQR_CS;IQR_NCS];
    TP=[TP_CS;TP_NCS];
    p100_200=[p100_200_CS;p100_200_NCS];
    p200_400=[p200_400_CS;p200_400_NCS];
    p400_800=[p400_800_CS;p400_800_NCS];
    spectrum_slope2=[spectrum_slope2_CS;spectrum_slope2_NCS];
    r_square2=[r_square2_CS;r_square2_NCS];
    
    meanPSD_label=[repmat(' meanPSD_CS', length(meanPSD_CS),1); repmat('meanPSD_NCS', length(meanPSD_NCS),1)];
    stdPSD_label=[repmat(' stdPSD_CS', length(stdPSD_CS),1); repmat('stdPSD_NCS', length(stdPSD_NCS),1)];
    medPSD_label=[repmat(' medPSD_CS', length(medPSD_CS),1); repmat('medPSD_NCS', length(medPSD_NCS),1)];
    bw_label=[repmat(' bw_CS', length(bw_CS),1); repmat('bw_NCS', length(bw_NCS),1)];
    p25_label=[repmat('p25_CS', length(p25_CS),1); repmat('p25NCS', length(p25_NCS),1)];
    p75_label=[repmat(' p75_CS', length(p75_CS),1); repmat('p75_NCS', length(p75_NCS),1)];
    IQR_label=[repmat(' IQR_CS', length(IQR_CS),1); repmat('IQR_NCS', length(IQR_NCS),1)];
    TP_label=[repmat(' TP_CS', length(TP_CS),1); repmat('TP_NCS', length(TP_NCS),1)];
    p100_200_label=[repmat(' p100_200_CS', length(p100_200_CS),1); repmat('p100_200_NCS', length(p100_200_NCS),1)];
    p200_400_label=[repmat(' p200_400_CS', length(p200_400_CS),1); repmat('p200_400_NCS', length(p200_400_NCS),1)];
    p400_800_label=[repmat(' p400_800_CS', length(p400_800_CS),1); repmat('p400_800_NCS', length(p400_800_NCS),1)];
    spectrum_slope2_label=[repmat(' spectrum_slope2_CS', length(spectrum_slope2_CS),1); repmat('spectrum_slope2_NCS', length(spectrum_slope2_NCS),1)];
    r_square2_label=[repmat(' r_square2_CS', length(r_square2_CS),1); repmat('r_square2_NCS', length(r_square2_NCS),1)];
    
    figure,
    subplot(w, h,1); boxplot(meanPSD, meanPSD_label,'Labels',{'meanPSD_CS', 'meanPSD_NCS'}); title('meanPSD');
    subplot(w, h,2); boxplot(stdPSD, stdPSD_label,'Labels',{'stdPSD_CS', 'stdPSD_NCS'}); title('stdPSD');
    subplot(w, h,3); boxplot(medPSD, medPSD_label,'Labels',{'medPSD_CS', 'meanPSD_NCS'}); title('medPSD');
    subplot(w, h,4); boxplot(bw,bw_label,'Labels',{'bw_CS', 'bw_NCS'}); title('bw');
    subplot(w, h,5); boxplot( p25,  p25_label,'Labels',{'p25_CS', 'p25_NCS'});title('p25');
    subplot(w, h,6); boxplot(p75, p75_label,'Labels',{'p75_CS', 'p75_NCS'});title('p75');
    subplot(w, h,7); boxplot(IQR, IQR_label,'Labels',{'IQR_CS', 'IQR_NCS'});title('IQR');
    subplot(w, h,8); boxplot(TP, TP_label,'Labels',{'TP_CS', 'TP_NCS'});title('TP');
    subplot(w, h,9); boxplot(p100_200, p100_200_label,'Labels',{'p100_200_CS', 'p100_200_NCS'});title('p100\_200');
    subplot(w, h,10); boxplot(p200_400, p200_400_label,'Labels',{'p200_400_CS', 'p200_400_NCS'});title('p200\_400');
    subplot(w, h,11); boxplot(p400_800, p400_800_label,'Labels',{'p400_800_CS', 'p400_800_NCS'});title('p400\_800');
    subplot(w, h,12); boxplot(spectrum_slope2, spectrum_slope2_label,'Labels',{'spectrum_slope2_CS', 'spectrum_slope2_NCS'});title('spectrum\_slope2');
    subplot(w, h,13); boxplot(r_square2, r_square2_label,'Labels',{'r_square2_CS', 'r_square2_NCS'});title('r\_square2');