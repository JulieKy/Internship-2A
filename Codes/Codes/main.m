%% Main script for spectral features analysis


%% INITIALISATION

clear all, close all, clc, dbstop if error;
addpath(genpath('..\..\')); % to have access to sample folder
init = 0; % optional, to generate a new Excel File
excelFile = 'NewFeaturesAnalysis_f'; % name of the Excel file to store features
path = pwd; % current path
pathExcel = strcat(path, '\'); % path of the Excel file

if ~(exist(pathExcel)) % test to create excel file or no
    disp('creation Excel folder')
    mkdir(pathExcel);
end

excelTemp = strcat([pathExcel excelFile], '.xls'); % add the .xls to have complete name
data_dir=[path,'\..\Data\Samples_Belle\'];

%% SPECTRAL FEATURES

%% Add the headers in the Excel file
if (init == 0)
    %     xlswrite([pathExcel excelFile], [{'file'}, {'ZCR'}, {'meanPSD'},{'stdPSD'},{'medPSD'},{'bw'},{'p25'},{'p75'},{'IQR'},{'TP'},{'p100_200'},{'p200_400'},{'p400_800'},{'SL'},{'R2'},{'nb_pks_MAF'}, {'f_higherPk_MAF'}, {'dif_higherPks_MAF'},{'nb_pks_GMM'}, {'f_higherPk_GMM'}, {'dif_higherPks_GMM'}, {'fi.a1'}, {'fi.a2'}, {'fi.a3'}, {'fi.a4'}, {'fi.b1'}, {'fi.b2'}, {'fi.b3'}, {'fi.b4'}, {'fi.c1'}, {'fi.c2'}, {'fi.c3'}, {'fi.c4'}, {'MFCC1'}, {'MFCC2'}, {'MFCC3'}, {'MFCC4'}, {'MFCC5'}, {'MFCC6'}], 'Features 1', 'A1'); % longest segment
    xlswrite([pathExcel excelFile], [{'file'}, {'ZCR'}], 'Temporal Features', 'A1'); % Sheet 1
    xlswrite([pathExcel excelFile], [{'file'}, {'meanPSD'},{'stdPSD'},{'medPSD'},{'bw'},{'p25'},{'p75'},{'IQR'},{'TP'},{'p100_200'},{'p200_400'},{'p400_800'},{'SL'},{'R2'}], 'Spectral Features 1', 'A1'); % Sheet 2
    xlswrite([pathExcel excelFile], [{'file'}, {'nb_pks_MAF'}, {'f_higherPk_MAF'}, {'dif_higherPks_MAF'},{'nb_pks_GMM'}, {'f_higherPk_GMM'}, {'dif_higherPks_GMM'}, {'fi.a1'}, {'fi.a2'}, {'fi.a3'}, {'fi.a4'}, {'fi.b1'}, {'fi.b2'}, {'fi.b3'}, {'fi.b4'}, {'fi.c1'}, {'fi.c2'}, {'fi.c3'}, {'fi.c4'}], 'Spectral Features 2', 'A1'); % Sheet 3
    xlswrite([pathExcel excelFile], [{'file'}, {'MFCC1'}, {'MFCC2'}, {'MFCC3'}, {'MFCC4'}, {'MFCC5'}, {'MFCC6'}, {'LPC2'}, {'LPC3'}, {'LPC4'}, {'LPC5'}, {'LPC6'}, {'LSF1'}, {'LSF2'}, {'LSF3'}, {'LSF4'}, {'LSF5'}, {'LSF6'} ], 'Coefficients', 'A1'); % Sheet 4
    init = 1;
end


%% Add the spectral features

%% Initialisation
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
for i = 1:lengthTot % loop to have all recording
    
    close all; % to close previous figures opened in computation of spectral features
    tempName=names_cell{i};
    % to follow which file is red
    disp('READ');
    disp(tempName);
    
    [x,Fs]= audioread([path,'\..\Data\Samples_Belle\',tempName]); % read current file
    
    %% resampling to 4000 Hz
    xs=resample(x,4000,Fs);
    fn=4000;
    
    %% Removing crying sections
    % -- Data
    observators=2;
    samples=37;
    end_sample=60; % End of the signal (hypotesis: length of the signal=60s)
    
    % -- Parameters
    window=3;
    overlap=25/100;
    
    % -- Function
    [labels, coef_KAPPA]=labelling(observators,samples, end_sample, window, overlap);
    
    %% filtering BP 100-1000Hz
    y = filterbp(xs,fn);
    
    %% Computation of features
    output_temporal_features = temporal_features(xs,fn, tempName); % Temporal features
    [output_spectral_features(i,:),periodogram_pks_features(i,:),pxx(i,:),f(i,:),foct(i,:),spower(i,:),I(i,:),S(i,:)] = spectral_features(y,fn); % See Fae's comment
    output_mean_mfcc = mfcc_coeffs(y, fn); % MFCCs coefficient
    [output_lpc, output_lsf] = lpc_lsf_coeff(y, fn); % LPC and LFC coefficient
    
    
    % -- Write on Excel file all the features
    
    % Sheet 1
    xlswrite([pathExcel excelFile], [i;output_temporal_features]', 'Temporal Features', ['A',num2str(i+1)]);
    xlswrite([pathExcel excelFile],{names_cell{i}}, 'Temporal Features',['A',num2str(i+1)]);
    
    % Sheet 2
    xlswrite([pathExcel excelFile], [i;output_spectral_features(i,:)']', 'Spectral Features 1', ['A',num2str(i+1)]);
    xlswrite([pathExcel excelFile],{names_cell{i}}, 'Spectral Features 1',['A',num2str(i+1)]);
    
    % Sheet 3
    xlswrite([pathExcel excelFile], [i;periodogram_pks_features(i,:)']', 'Spectral Features 2', ['A',num2str(i+1)]);
    xlswrite([pathExcel excelFile],{names_cell{i}}, 'Spectral Features 2',['A',num2str(i+1)]);
    
    % Sheet 4
    xlswrite([pathExcel excelFile], [i;output_mean_mfcc'; output_lpc'; output_lsf']', 'Coefficients', ['A',num2str(i+1)]);
    xlswrite([pathExcel excelFile],{names_cell{i}}, 'Coefficients',['A',num2str(i+1)]);
    
    %     spectrogram(y, 'yaxis')
end

%% Plot Average figures


%% FFT Representation

% Median with inter-quartile range
figure(),
plot(mean(f(find(f(:,2) ~= 0),:)) , median(pxx(find(f(:,2) ~= 0),:))); % average fft
hold on
plot(mean(f(find(f(:,2) ~= 0),:)) , prctile(pxx(find(f(:,2) ~= 0),:),25), 'LineStyle', '--', 'Color', 'r'); % std
plot(mean(f(find(f(:,2) ~= 0),:)) , prctile(pxx(find(f(:,2) ~= 0),:),75), 'LineStyle', '--', 'Color', 'r'); % std
xlabel('Frequency (in Hz)');
ylabel('Power Spectrum');
title('Average Power Spectrum (All groups)');
legend('Average power spectrum', 'Interquartile range');
hold off

% Mean plot
figure(),
plot(mean(f(find(f(:,2) ~= 0),:)) , mean(pxx(find(f(:,2) ~= 0),:)));
xlabel('Frequency (in Hz)');
ylabel('Power Spectrum');
title('Average Power Spectrum (all)');
legend('Average Power Spectrum');
