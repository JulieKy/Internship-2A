% Spectrogram of NCS and CS

%% INITIALISATION PATH
clear all, close all, clc, dbstop if error;
addpath(genpath('..\..\')); % to have access to sample folder
init = 0; % optional, to generate a new Excel File
excelFile = 'NewFeaturesAnalysis_f'; % name of the Excel file to store features
path = pwd; % current path
data_dir=[path,'\..\..\Data\Samples_Belle\'];

%% READING SIGNAL 15
tempName='15.mp3';
signal_n=15;
str=sprintf(' -- READ: %s --\n', tempName);
disp(str)
[y,Fs]= audioread([path,'\..\..\Data\Samples_Belle\',tempName]); % read current file

%% PREPROCESSING
% -- Resampling to 4000 Hz
xs=resample(y,4000,Fs);
fn=4000;

% -- Shorten to 60s
time_sample=60;
xss=xs(1:time_sample*fn,1);

% -- Parameters
N = length(xss);
time_axis = (1:N)/fn;

%% SPECTROGRAM
% -- Parameters A METTRE DANS APPEL FONCTION
wind_time=1; % Window of 1 second
overlap=0.25; % 25% overlap
start_time=0;
end_time=15;

% -- Parameters in samples
wind_sample= wind_time*fn;
overlap_sample=wind_sample*overlap; % 25% overlap
start_sample=start_time*fn+1;
end_sample=end_time*fn+1;

% -- Part of the signal wanted
y= xss(start_sample:end_sample); % First 15 second of signal 15.mp3
time_axis_y=time_axis(start_sample:end_sample);

figure,
subplot(2,1,1);
plot(time_axis_y, y, 'Color', [0 0.6 0]);
str=sprintf('Signal %d', signal_n); 
legend(str)
xlabel('Time (s)')
ylabel('Amplitude')
title('Time representation of CS and NCS')

% -- Display
subplot(2,1,2);
spectrogram(y,wind_sample,overlap_sample, 128, fn, 'yaxis');
str=sprintf('Spectrogram of the fisrt %d seconds of %d.mp3', end_time-start_time, signal_n); 
title(str)
suptitle('Vizualisation of CS frequency changes');




%% LABELLING

%% -- Initialisation
% -- Data
observators=3;
samples=37;
end_sample=60; % End of the signal (hypotesis: length of the signal=60s)

% -- Parameters for labelling
window=1;
overlap=0;

[label_final, coef_KAPPA]=labelling(observators,samples, end_sample, window, overlap);


signal_n=13;
display_NCS_CS_annotations(signal_n,label_final, window, overlap)

signal_n=15;
display_NCS_CS_annotations(signal_n,label_final, window, overlap, start_time, end_time)















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
% end
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





