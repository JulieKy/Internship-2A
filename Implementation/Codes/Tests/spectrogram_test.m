% Spectrogram of NCS and CS

%% INITIALISATION PATH
clear all, close all, clc, dbstop if error;
addpath(genpath('..\..\')); % to have access to sample folder
init = 0; % optional, to generate a new Excel File
excelFile = 'NewFeaturesAnalysis_f'; % name of the Excel file to store features
path = pwd; % current path
data_dir=[path,'\..\..\Data\Samples_Belle\'];
tempName='15.mp3';

signal_n=15;

str=sprintf(' -- READ: %s --\n', tempName);
disp(str)

[x,Fs]= audioread([path,'\..\..\Data\Samples_Belle\',tempName]); % read current file

%% -- Resampling to 4000 Hz
xs=resample(x,4000,Fs);
fn=4000;

%% -- Shorten to 60s
time_sample=60;
xss=xs(1:time_sample*fn,1);

%% INITIALISATION
% -- Data
observators=3;
samples=37;
end_sample=60; % End of the signal (hypotesis: length of the signal=60s)

% -- Parameters for labelling
window=1;
overlap=0;

%% LABELLING
[label_final, coef_KAPPA]=labelling(observators,samples, end_sample, window, overlap);

%% -- Finding NCS/CS on the signal
% -- Variables
N = length(xss);
time_axis = (1:N)/fn;

% % -- Finding the location of 'CS'
% flag_section=1; %CS
% locs=find(label_final(signal_n,:)==flag_section); % Locations of NCS/CS
% 
% if isempty(locs)==0 % There are NCS/CS on the signal
%     
%     % -- Start time of the labels (for each window)
%     start_time=locs*(window-window*overlap);
%     
%     % -- Start sample of the labels (for each window)
%     start_sample=start_time*fn;
%     label_duration=window*fn; % Number of samples in a window
%     
%     figure,
%     % for n_section=1:length(locs)
%     n_section=1
%     [xss_section,time_axis_section] = label2signal(xss, n_section, start_sample, label_duration, time_axis);
%     %         subplot(1,length(locs),n_section)
%     spectrogram(xss_section,  'yaxis')
%     str=sprintf('Spectrogram CS %d of 15.mp3',n_section);
%     title(str);
%     figure,
%     plot(time_axis, xss); hold on
%     plot(time_axis_section,xss_section); hold off
%     str=sprintf('CS %d',n_section);
%     title(str);
% end
% % end
% 
% % -- Finding the location of 'NCS'
% flag_section=0; %NCS
% locs=find(label_final(signal_n,:)==flag_section); % Locations of NCS/CS
% 
% if isempty(locs)==0 % There are NCS/CS on the signal
%     
%     % -- Start time of the labels (for each window)
%     start_time=locs*(window-window*overlap);
%     
%     % -- Start sample of the labels (for each window)
%     start_sample=start_time*fn;
%     label_duration=window*fn; % Number of samples in a window
%     
%     figure,
%     %     for n_section=1:length(locs)
%     [xss_section,time_axis_section] = label2signal(xss, n_section, start_sample, label_duration, time_axis);
%     %         subplot(1,length(locs),n_section)
%     spectrogram(xss_section,  'yaxis')
%     str=sprintf('Spectrogram NCS %d of 15.mp3',n_section);
%     title(str);
%     figure,
%     plot(time_axis, xss); hold on
%     plot(time_axis_section,xss_section); hold off
%     str=sprintf('NCS %d',n_section);
%     title(str);
%     %     end
% end
% 
% 
% figure,
% spectrogram(xss,  'yaxis')
% title('Spectrogram of the entire signal 15.mp3')


display_NCS_CS_annotations(signal_n,label_final, window, overlap)

