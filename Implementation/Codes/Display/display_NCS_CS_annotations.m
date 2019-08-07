function [] = display_NCS_CS_annotations(signal_n,label_final, window, overlap)
%NCS_CS_LABEL: Display the annotated labels CS and NCS of a signal

%% INPUTS AND OUTPUTS
%  -- Inputs --
% signal_n: Number of the signal
% label_final: Annotated labels of the signal bank


%% Reading the signal
signal_name=sprintf('%d.mp3',signal_n);
path=pwd;
[x1,Fs]= audioread([path,'\..\Data\Samples_Belle\',signal_name]); % read current file

%% Resampling to 4000 Hz
xs1=resample(x1,4000,Fs);
fn=4000;

%% Shorten the signals to 60s
time_sample=60;
xss1=xs1(1:time_sample*fn,1);

%% Parameters
N = length(xss1);
time_axis = (1:N)/fn;

NCS_color=[0.4 1 0.4];
CS_color=[1 0.2 0.3];

%% Display
figure,

% -- Temporal representation of the signal
subplot(2,1,1);
plot(time_axis, xss1);
title('Temporal Representation of Signal 22')
xlabel('Time [s]');
ylabel('Amplitude');

% -- Temporal representation of the signal with NCS and CS
subplot(2,1,2)
plot(time_axis, xss1, 'Color', CS_color); hold on % The entire signal %

% Finding the location of 'CS' and 'NCS'
CS_locs=find(label_final(signal_n,:)==1);
NCS_locs=find(label_final(signal_n,:)==0);

% Start time of the labels (for each window)
time_CS_start=CS_locs*(window-overlap);
time_NCS_start=NCS_locs*(window-overlap);

% Start sample of the labels (for each window)
sample_CS_start=time_CS_start*fn;
sample_NCS_start=time_NCS_start*fn;
label_duration=window*fn; % Number of samples in a window

% For each NCS section
for n_section=1:length(NCS_locs)
    NCS_start=sample_NCS_start(n_section);
    NCS_end=NCS_start+label_duration;
    
    % Edge effects
    if NCS_end>length(xss1) % For the last section, if the window and overlap or not adjusted to xss length
        NCS_end=length(xss1);
    end
    
    % Section on the signal
    NCS_section=xss1(NCS_start:NCS_end);
    time_axis_section=time_axis(NCS_start:NCS_end);
    plot(time_axis_section, NCS_section, 'Color', NCS_color);
end

% For each CS section
for n_section=1:length(CS_locs)
    CS_start=sample_CS_start(n_section);
    CS_end=CS_start+label_duration;
    
    % Edge effects
    if CS_end>length(xss1) % For the last section, if the window and overlap or not adjusted to xss length
        CS_end=length(xss1);
    end
    
    % Section on the signal
    CS_section=xss1(CS_start:CS_end);
    time_axis_section=time_axis(CS_start:CS_end);
    plot(time_axis_section, CS_section, 'Color', CS_color);
end

hold off
legend('NCS', 'CS')
title('CS and NCS after Annotations on Signal 22')
xlabel('Time [s]');
ylabel('Amplitude');

end

