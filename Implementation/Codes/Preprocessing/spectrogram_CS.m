function [s, f, t] = spectrogram_CS(signal_n, wind_time_spec, overlap_spec, start_time, end_time)
%display_spectrogram_CS: Display the spectrogram and the time
%representation with CS and NCS of a part of a signal AND return its
%spectrogram

%% INPUTS AND OUTPUTS
%  -- Inputs --
% signal_n: Signal number
% wind_time_spec: Duration of the spectrogram window
% overlap_spec: Spectrogram overlap
% start_time: Beginning of the part wanted (s)
% end_time: End of the part wanted (s)
%  -- Outputs --
% [s, f, t]: Spectrogram of the signal part wanted

%% INITIALISATION PATH
addpath(genpath('..\..\')); % to have access to sample folder
init = 0; % optional, to generate a new Excel File
excelFile = 'NewFeaturesAnalysis_f'; % name of the Excel file to store features
path = pwd; % current path
data_dir=[path,'\..\Data\Samples_Belle\'];

%% READING SIGNAL 15
tempName=sprintf('%d.mp3', signal_n);
str=sprintf(' -- READ: %s --\n', tempName);
disp(str)
[y,Fs]= audioread([path,'\..\Data\Samples_Belle\',tempName]); % read current file

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

%% LABELLING
% -- Data
observators=3;
samples=37;
end_sample=60; % End of the signal (hypotesis: length of the signal=60s)

% -- Parameters for labelling
window=1;
overlap=0;

% -- Labels
[label_final, coef_KAPPA]=labelling(observators,samples, end_sample, window, overlap);


%% PARAMETERS

% -- Parameters in samples
wind_sample_spec= wind_time_spec*fn;
overlap_sample_spec=wind_sample_spec*overlap_spec; % 25% overlap

% -- Part of the signal wanted
start_sample=start_time*fn+1;
end_sample=end_time*fn+1;
y= xss(start_sample:end_sample); % First 15 second of signal 15.mp3
time_axis_y=time_axis(start_sample:end_sample);


%% CS LOCATIONS
% -- Finding the location of 'CS'
flag_section=1; %CS
CS_locs=find(label_final(signal_n,:)==flag_section); % Locations of CS

if isempty(CS_locs)==0 % There are CS on the signal
    
    % -- Start time of the labels (for each window)
    start_time_CS=CS_locs*(window-window*overlap)-1;
    
    % -- Start sample of the labels (for each window)
    sample_CS_start=start_time_CS*fn+1;
    label_duration=window*fn; % Number of samples in a window
end

[s, f, t]= spectrogram(y,wind_sample_spec,overlap_sample_spec, 128, fn, 'yaxis');

%% DISPLAY

figure,

% -- Display NCS/CS
subplot(2,1,1);
p1=plot(time_axis_y, y, 'Color', [0 0.6 0]);hold on

for n_section=1:length(CS_locs)
    if CS_locs(n_section)>=start_time && CS_locs(n_section)<=end_time
        [CS_section,time_axis_section] = label2signal(y, n_section, sample_CS_start, label_duration, time_axis);
        p2=plot(time_axis_section, CS_section, 'Color', [0.8 0 0]);
    end
end
hold off;
legend([p1 p2],'NCS', 'CS')
xlabel('Time (s)')
ylabel('Amplitude')
str=sprintf('Time Representation of the fisrt %ds of %d.mp3', end_time-start_time, signal_n);
title(str)

% -- Display spectrogram
subplot(2,1,2);
spectrogram(y,wind_sample_spec,overlap_sample_spec, 128, fn, 'yaxis');
str=sprintf('Spectrogram of the fisrt %ds of %d.mp3', end_time-start_time, signal_n);
title(str)
str2=sprintf('Vizualisation of CS Frequency Changes in Signal %d.mp3',signal_n);
suptitle(str2);

end

