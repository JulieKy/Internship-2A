function [X_ncs_final, label_training] = crying_removing(X, fn, threshold, band, window_training, overlap_label, window_annotated)
%CRYING_REMOVING:  Remove the crying sections thanks to a threshold on powerband

%% INPUTS AND OUTPUTS
%  -- Inputs --
% X: Matrix with signals
% fn: Sampling frequency
% threshold: Power band threshold found to distinguish CS from NCS
% band: Frequency band of interest for doing the power ratio and distinguish CS from NCS
% window_training: Window used for training
% window_annotated: Window used for annotated
% -- Outputs --
% X_ncs_final: Matrix with non-crying signals
% label_training: Label found by training (useful in display_CS_NCS_final.m, figure 3)

X_ncs=zeros(size(X,1), size(X,2));

for signal_n=1:size(X,1)
    xss=X(signal_n, :);
    
    %% INITIALISATION
    N=length(xss); % Signal length
    duration_sample=window_training*fn; % Duration of the section in sample
    nb_section=floor(N/duration_sample); % Number of sections
    n_section=1:nb_section; % Section number (ID)
    start_time=n_section*(window_training-window_training*overlap_label)-1; % Start time of each section
    start_sample=start_time*fn+1; % Start sample of each section
    
    %% POWER BAND OF EACH SECTION
    xss_section=reshape(xss,  duration_sample, [length(xss)/duration_sample]); % Each section in a column
    power_band=bandpower(xss_section, fn, band); % Power band of each section
    
    %% CS and NCS
    CS=power_band>threshold;
    NCS=1-CS;   
    
    %% CS REMOVING
    xss_section_NCS=xss_section.*NCS; % Each CS column is filled with 0
    xss_NCS=reshape(xss_section_NCS, [1 size(xss_section_NCS,1)*size(xss_section_NCS,2)]); % Reshape in a single line
    
    xsc=xss_NCS(xss_NCS~=0); % Signal without CS
    
    %% STORAGE
    X_ncs(signal_n, 1:length(xsc))=xsc; % If CS were removed, X_ncs contains 0.
    
    %% MINIMUM LENGTH
    length_xsc(signal_n)=length(xsc);
    
    %% LABEL TRAINING
    label_training(signal_n, :)=repelem(CS, round(window_training/window_annotated)); % Useful in 'display_CS_NCS_final.m'
    
end

%% SHORTEN SAMPLES WITH MINIMUM LENGTH
% Find the minimum length
min_length=min(length_xsc);
part1=floor(min_length/2);
part2=min_length-part1;

% Statistical study of lengths
length_xsc_time=length_xsc*60/length(xss);
mean_length=mean(length_xsc_time); 
median_length=median(length_xsc_time);
p25_length=prctile(length_xsc_time,25);
p75_length=prctile(length_xsc_time,75);


% Shorten for each signal
X_ncs_final=zeros(size(X,1), min_length);
for signal_n=1:size(X_ncs,1)
    xsc1=X_ncs(signal_n, 1:length_xsc(signal_n)); % Signal with only NCS
    if length(xsc1)>min_length
        mid=floor(length(xsc1)/2);
        xsc_shorten=xsc1(mid-part1:mid+part2-1);
    else 
        xsc_shorten=xsc1;
    end
    X_ncs_final(signal_n, :)=xsc_shorten;
end

end
