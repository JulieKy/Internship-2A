%% Main script for spectral features analysis


%% INITIALISATION

clear all, close all, clc;
addpath(genpath('../')); % to have access to sample folder
init = 0; % optional, to generate a new Excel File
excelFile = 'NewFeaturesAnalysis_f'; % name of the Excel file to store features
path = pwd; % current path
pathExcel = strcat(path, '\'); % path of the Excel file

if ~(exist(pathExcel)) % test to create excel file or no
    disp('creation Excel folder')
    mkdir(pathExcel);
end

excelTemp = strcat([pathExcel excelFile], '.xls'); % add the .xls to have complete name
data_dir=[path,'\Samples\'];

%% SPECTRAL FEATURES

%% Add the headers in the Excel file
if (init == 0)
    xlswrite([pathExcel excelFile], [{'file'}, {'meanPSD'},{'stdPSD'},{'medPSD'},{'bw'},{'p25'},{'p75'},{'IQR'},{'TP'},{'p100_200'},{'p200_400'},{'p400_800'},{'SL'},{'R2'}, {'MFCC1'}, {'MFCC2'}, {'MFCC3'}, {'MFCC4'}, {'MFCC5'}, {'MFCC6'}], 'Features 1', 'A1'); % longest segment
      init = 1;
end


%% Add the spectral features

%% Initialisation
dinfo = dir(data_dir);
names_cell1 = {dinfo.name};
    % choose a valid file name
    j=0;
 for i=1:size(names_cell1,2)
       if length(names_cell1{i})>4
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
    
    [x,Fs]= audioread([path,'\Samples\',tempName]); % read current file
    
    %% resampling to 4000 Hz
      xs=resample(x,4000,Fs);
    fn=4000;
    %% filtering BP 100-1000Hz
    y = filterbp(xs,fn);
    %% Computation of features
    output_mean_mfcc = mfcc_coeffs(y, fn); % MFCCs coefficient 
        [output_spectral_features(i,:),pxx(i,:),f(i,:),foct(i,:),spower(i,:),I(i,:),S(i,:)] = spectral_features(y,fn); % See Fae's comment         
        % write on Excel file all the features
        xlswrite([pathExcel excelFile], [i;output_spectral_features(i,:)' ; output_mean_mfcc']', 'Features 1', ['A',num2str(i+1)]); % Sheet 1
    xlswrite([pathExcel excelFile],{names_cell{i}}, 'Features 1',['A',num2str(i+1)]); 
    
    
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
