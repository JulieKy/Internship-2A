function [pxx_mean, band_mean, PR_mean, f, p25_signal, p75_signal]=power_ratio_band(xss, signal_n, fn, window, overlap, label_final, pass_band, band_width, flag_section)
%POWER RATIO:  Calculate the power ratio of NCS/CS (depending on flag_section) of a signal

%% INPUTS AND OUTPUTS
%  -- Inputs --
%
% -- Outputs --


%% INITIALISATION

% -- Variables
N = length(xss);
time_axis = (1:N)/fn;

% -- Finding the location of 'CS' or 'NCS'
locs=find(label_final(signal_n,:)==flag_section); % Locations of NCS/CS

if isempty(locs)==1 % There isn't NCS/CS on the signal
    band_mean=0;
    pxx_mean=0;
    PR_mean=0;
    f=0;
    p25_signal=0;
    p75_signal=0;
    
else % There are NCS/CS on the signal
    
    % -- Start time of the labels (for each window)
    start_time=locs*(window-window*overlap)-1;
    
    % -- Start sample of the labels (for each window)
    start_sample=start_time*fn+1;
    label_duration=window*fn; % Number of samples in a window
    
    
    %% PERIODOGRAM
    
    pxx_signal=[]; % Stores the periodogram of every NCS/CS of a signal
    p25_signal=[]; % Stores the first interquartile of every NCS/CS of a signal
    p75_signal=[]; % Stores the third interquartile of every NCS/CS of a signal
    
    % -- For each section (NCS/CS)
    for n_section=1:length(locs)
        
        % Section on the signal (NCS/CS)
        [section,time_axis_section] = label2signal(xss, n_section, start_sample, label_duration, time_axis);
        
        % Welch Periodogram of the section
        [pxx,f] = Welch_periodogram(section, fn, pass_band); % Welch Periodogram of the section
        pxx_signal=[pxx, pxx_signal];
 
       
        % Interquartile in frequency of each section
        p25=percentilefreq(pxx((f~=0)),f(f~=0),25);
        p75=percentilefreq(pxx((f~=0)),f(f~=0),75);
        p25_signal=[p25, p25_signal];
        p75_signal=[p75, p75_signal];
       
        % Parameters
        f_interval=length(f)*band_width/(pass_band(end)-pass_band(1)); % Interval of frequencies in the vector f
        end_band=length(f);
        band_mean_section=[]; % Means of the periodogram in the frequency bands
        
        % For each frequency band
        for n_band = 1 : length(f)/f_interval
            start_band=end_band-(floor(f_interval)); % Beginning of the band
            pxx_band=pxx(start_band:end_band); % Periodogram in the band
            band_mean_section=[ band_mean_section, mean(pxx_band)]; % Mean of this periodogram
            end_band=start_band; % End of the band
        end
        
        band_mean_signal(n_section,:)=band_mean_section; % Periodogram means in each frequency bands, for each section NCS/CS of a signal. Matrix shaped like this: (section, frequency bands)
        PR=band_mean_signal/sum(pxx); % Power ratio in each frequency bands, for each section NCS/CS of a signal
        
    end
    
    
    %% OUTPUTS 
    if length(locs)>1 % More than one section
        band_mean=mean(band_mean_signal); % For a signal, mean of the means of frequency bands periodogram, for all NCS/CS sections
        pxx_mean=mean(pxx_signal'); % For a signal, mean of the periodograms of all NCS/CS sections
        PR_mean=mean(PR); % For a signal, mean of power ratios of all NCS/CS sections
    else
        band_mean=band_mean_signal; % For a signal, mean of the means of frequency bands periodogram, for all NCS/CS sections
        pxx_mean=pxx_signal'; % For a signal, mean of the periodograms of all NCS/CS sections
        PR_mean=PR; % For a signal, mean of power ratios of all NCS/CS sections
    end
end


