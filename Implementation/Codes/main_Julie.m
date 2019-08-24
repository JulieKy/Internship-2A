% MAIN_JULIE: Main script of the project (preprocessing + spectral analysis)

%% ------ INITIALISATION ------

clear all, close all, clc, dbstop if error;

%% -- CS and NCS color
NCS_color=[0 0.6 0];
CS_color=[0.8 0 0];

%% -- Signal parameters
time_sample=60; % Signal duration
fn=4000; % Sampling frequency

%% -- Path
addpath(genpath('..\..\')); % Access to sample folder
path = pwd; % Current path
data_dir=[path,'\..\Data\Samples_Belle\'];
pathExcel = strcat(path, '\'); % Path of the Excel spectral features file
pathExcelPreprocessing = strcat(path, '\'); % Path of the Excel preprocessing file

% %% -- Excel file Spectral Features
% init = 0; % To generate a new Excel files if necessary
% excelFileSpectralFeatures = 'Spectral_Features'; % Name of the Excel file to store features
% 
% if ~(exist(pathExcel)) % test to create excel file or no
%     disp('creation Excel Spectral Features folder')
%     mkdir(pathExcel);
% end
% 
% excelTemp = strcat([pathExcel excelFileSpectralFeatures], '.xls'); % add the .xls to have complete name
% 
% % Add the headers in the Excel file
% if (init == 0)
%     xlswrite([pathExcel excelFileSpectralFeatures], [{'file'}, {'ZCR'}, {'meanPSD'},{'stdPSD'},{'medPSD'},{'bw'},{'p25'},{'p75'},{'IQR'},{'TP'},{'p100_200'},{'p200_400'},{'p400_800'},{'SL'},{'R2'},{'nb_pks_MAF'}, {'f_higherPk_MAF'}, {'dif_higherPks_MAF'},{'nb_pks_GMM'}, {'f_higherPk_GMM'}, {'dif_higherPks_GMM'}, {'fi.a1'}, {'fi.a2'}, {'fi.a3'}, {'fi.a4'}, {'fi.b1'}, {'fi.b2'}, {'fi.b3'}, {'fi.b4'}, {'fi.c1'}, {'fi.c2'}, {'fi.c3'}, {'fi.c4'}, {'MFCC1'}, {'MFCC2'}, {'MFCC3'}, {'MFCC4'}, {'MFCC5'}, {'MFCC6'}], 'Features 1', 'A1'); % longest segment
%     xlswrite([pathExcel excelFileSpectralFeatures], [{'file'}, {'ZCR'}], 'Temporal Features', 'A1'); % Sheet 1
%     xlswrite([pathExcel excelFileSpectralFeatures], [{'file'}, {'meanPSD'},{'stdPSD'},{'medPSD'},{'bw'},{'p25'},{'p75'},{'IQR'},{'TP'},{'p100_200'},{'p200_400'},{'p400_800'},{'SL'},{'R2'}], 'Spectral Features 1', 'A1'); % Sheet 2
%     xlswrite([pathExcel excelFileSpectralFeatures], [{'file'}, {'nb_pks_MAF'}, {'f_higherPk_MAF'}, {'dif_higherPks_MAF'},{'nb_pks_GMM'}, {'f_higherPk_GMM'}, {'dif_higherPks_GMM'}, {'fi.a1'}, {'fi.a2'}, {'fi.a3'}, {'fi.a4'}, {'fi.b1'}, {'fi.b2'}, {'fi.b3'}, {'fi.b4'}, {'fi.c1'}, {'fi.c2'}, {'fi.c3'}, {'fi.c4'}], 'Spectral Features 2', 'A1'); % Sheet 3
%     xlswrite([pathExcel excelFileSpectralFeatures], [{'file'}, {'MFCC1'}, {'MFCC2'}, {'MFCC3'}, {'MFCC4'}, {'MFCC5'}, {'MFCC6'}, {'LPC2'}, {'LPC3'}, {'LPC4'}, {'LPC5'}, {'LPC6'}, {'LSF1'}, {'LSF2'}, {'LSF3'}, {'LSF4'}, {'LSF5'}, {'LSF6'} ], 'Coefficients', 'A1'); % Sheet 4
%     init = 1;
% end
% 
% %% -- Excel file Preprocessing
 init_learning=0;
excelFileSpectralFeaturesPreprocessing = 'CS_Features'; % Name of the Excel file to store features

if ~(exist(pathExcelPreprocessing)) % test to create excel file or no
    disp('Creation Excel Preprocessing folder')
    mkdir(pathExcelPreprocessing);
end

excelTempPreprocessing = strcat([pathExcelPreprocessing excelFileSpectralFeaturesPreprocessing], '.xls'); % add the .xls to have complete name

% Add the headers in the Excel file
if (init_learning == 0)
    xlswrite([pathExcelPreprocessing excelFileSpectralFeaturesPreprocessing], [{'Threshold'},  {'p25'}, {'p75'}, {'Window_label'}, {'Overlap_label'}], 'Learning CS Features', 'A1');
end



%% -- Samples' names initialisation
dinfo = dir(data_dir);
names_cell1 = {dinfo.name};

% Choose a valid file name
j=0;
for i=1:size(names_cell1,2)
    if length(names_cell1{i})>2
        j=j+1;
        names_cell{j}=names_cell1{i};
    end
end
lengthTot=j;

%% -- STORAGE OF SAMPLES
X=zeros(lengthTot, time_sample*fn); % All the signals stored in this matrix

for i = 1:lengthTot % loop to have all recording
    
    close all; % Close previous figures opened in computation
    
    % Name of the sample
    tempName=names_cell{i};
    disp('READ - main_Julie.m');
    disp(tempName);
    
    % Get the number of the recording by removing the '.mp3'
    strMP3 = sprintf('%s',tempName);
    ind=strfind(strMP3,'.');
    signal_n = str2num(strMP3(1:ind-1));
    
    % Reading the sample
    [x,Fs]= audioread([path,'\..\Data\Samples_Belle\',tempName]); % read current file
    
    % Resampling to 4000 Hz
    xs=resample(x,fn,Fs);
    
    % Shorten the signals to 60s
    xss=xs(1:time_sample*fn,1);
    
    % Fill X with all samples
    X(signal_n, :)=xss;
end

%% ------  PREPROCESSING ------

%% -- Removing crying sections (CS)
overlap_label=0;

% -- Learning where are the CS (if not already done)
if (init_learning == 0) % Need to be done one time (data stored on an Excel file)
    [threshold,  band, label_annotated, window_annotated, window_training]= crying_learning(X, fn, CS_color, NCS_color);
    xlswrite([pathExcelPreprocessing excelFileSpectralFeaturesPreprocessing], [threshold ; band(1); band(end); window_annotated; window_training]', 'Learning CS Features', 'A2');
    xlswrite([pathExcelPreprocessing excelFileSpectralFeaturesPreprocessing], [label_annotated], 'Annotated Labels', 'A1');
    init_learning=1;
end

% -- Removing the CS
% Reading data from the Excel file
outputs_ExcelProcessing=xlsread([pathExcelPreprocessing excelFileSpectralFeaturesPreprocessing],'Learning CS Features','A2:E2');
threshold=outputs_ExcelProcessing(1); band(1)=outputs_ExcelProcessing(2); band(2)=outputs_ExcelProcessing(3); window_annotated=outputs_ExcelProcessing(4); window_training=outputs_ExcelProcessing(5);
label_annotated=xlsread([pathExcelPreprocessing excelFileSpectralFeaturesPreprocessing], 'Annotated Labels');

% Removing the data
[X_ncs, label_learning]=crying_removing(X, fn, threshold, band, window_training, overlap_label, window_annotated);


%% -- Display NCS and CS
signal_n=1;
xss=X(signal_n,:);
xsc=X_ncs(signal_n,:);

% Display xss, annotated labels, learnt labels and xsc
display_CS_NCS_final(xss, xsc, fn, signal_n, label_annotated, window_annotated, overlap_label, label_learning, NCS_color, CS_color);


%% -- Filtering BP 100-1000Hz
y = filterbp(xsc,fn);

%% ------  SPECTRAL FEATURES ------

%% -- Computation of features
output_temporal_features = temporal_features(xs,fn, tempName); % Temporal features
[output_spectral_features(i,:),periodogram_pks_features(i,:),pxx(i,:),f(i,:),foct(i,:),spower(i,:),I(i,:),S(i,:)] = spectral_features(y,fn); % See Fae's comment
output_mean_mfcc = mfcc_coeffs(y, fn); % MFCCs coefficient
[output_lpc, output_lsf] = lpc_lsf_coeff(y, fn); % LPC and LFC coefficient


%% -- Write on Excel file all the features

% Sheet 1
xlswrite([pathExcel excelFileSpectralFeatures], [i;output_temporal_features]', 'Temporal Features', ['A',num2str(i+1)]);
xlswrite([pathExcel excelFileSpectralFeatures],{names_cell{i}}, 'Temporal Features',['A',num2str(i+1)]);

% Sheet 2
xlswrite([pathExcel excelFileSpectralFeatures], [i;output_spectral_features(i,:)']', 'Spectral Features 1', ['A',num2str(i+1)]);
xlswrite([pathExcel excelFileSpectralFeatures],{names_cell{i}}, 'Spectral Features 1',['A',num2str(i+1)]);

% Sheet 3
xlswrite([pathExcel excelFileSpectralFeatures], [i;periodogram_pks_features(i,:)']', 'Spectral Features 2', ['A',num2str(i+1)]);
xlswrite([pathExcel excelFileSpectralFeatures],{names_cell{i}}, 'Spectral Features 2',['A',num2str(i+1)]);

% Sheet 4
xlswrite([pathExcel excelFileSpectralFeatures], [i;output_mean_mfcc'; output_lpc'; output_lsf']', 'Coefficients', ['A',num2str(i+1)]);
xlswrite([pathExcel excelFileSpectralFeatures],{names_cell{i}}, 'Coefficients',['A',num2str(i+1)]);


%% ****** DISPLAY ******

%% -- Display NCS and CS
% CS and NCS color
NCS_color=[0 0.6 0];
CS_color=[0.8 0 0];

% Display xss, annotated labels, learnt labels and xsc
display_CS_NCS_final(xss, xsc, fn, signal_n, label_annotated, window_label, overlap_label, label_learning_xss, NCS_color, CS_color);


%% -- FFT Representation

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
