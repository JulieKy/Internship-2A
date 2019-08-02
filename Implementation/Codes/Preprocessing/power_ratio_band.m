function [pxx_mean, band_mean, PR_mean, f]=power_ratio_band(xss, fn, window, overlap, label_final, pass_band, band_width, flag_section)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

i=1; % NUMERO DE SAMPLE A DONNER!! OU VOIR, ON l'a dans l'qutre???

%% INITIALISATION

% -- Variables
N = length(xss);
time_axis = (1:N)/fn;

% -- Finding the location of 'CS' or 'NCS'
locs=find(label_final(i,:)==flag_section); % VOIR LE i

% -- Start time of the labels (for each window)
start_time=locs*(window-window*overlap);

% -- Start sample of the labels (for each window)
start_sample=start_time*fn;
label_duration=window*fn; % Number of samples in a window


%% PERIODOGRAM

pxx_signal=[]; % Stores the periodogram of every NCS/CS of a signal

% -- For each section (NCS/CS)
for n_section=1:length(locs)
    
    % Beginning and end samples of a section
    start_section=start_sample(n_section);
    end_section=start_section+label_duration;
    
    % Edge effects
    if end_section>length(xss) % For the last section, if the window and overlap or not adjusted to xss length
        end_section=length(xss);
    end
    
    % Section on the signal (NCS/CS)
    section=xss(start_section:end_section);
    time_axis_section=time_axis(start_section:end_section);
     
    % Welch Periodogram of the section
    [pxx,f] = Welch_periodogram(section, fn, pass_band); % Welch Periodogram of the section
    pxx_signal=[pxx, pxx_signal]; 
    
    % Parameters
    f_interval=length(f)*band_width/(pass_band(end)-pass_band(1)); % Interval of frequencies in the vector f
    end_band=length(f);
    band_mean_section=[]; % Means of the periodogram in the frequency bands
    
    % For each frequency band
    for n_band = 1 : length(f)/f_interval 
        start_band=end_band-(floor(f_interval)); % Beginning of the band
        pxx_band=pxx(start_band:end_band); % Periodogram in the band
        band_mean_section=[mean(pxx_band), band_mean_section]; % Mean of this periodogram
        end_band=start_band; % End of the band     
    end
    
    band_mean_signal(n_section,:)=band_mean_section; % Periodogram means in each frequency bands, for each section NCS/CS of a signal. Matrix shaped like this: (section, frequency bands)
    PR=band_mean_signal/sum(pxx); % Power ratio in each frequency bands, for each section NCS/CS of a signal
end

%% OUTPUTS
band_mean=mean(band_mean_signal); % For a signal, mean of the means of frequency bands periodogram, for all NCS/CS sections
pxx_mean=mean(pxx_signal'); % For a signal, mean of the periodograms of all NCS/CS sections
PR_mean=mean(PR); % For a signal, mean of power ratios of all NCS/CS sections


end

