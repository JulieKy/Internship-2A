function [xsc, CS] = crying_removing(xss, fn, threshold, band, window_label, overlap_label)
%CRYING_REMOVING:  Remove the crying sections thanks to a threshold on powerband

%% INPUTS AND OUTPUTS
%  -- Inputs --
% xs: signal
% -- Outputs --
% xsc: signal without crying sections

%% INITIALISATION
N=length(xss); % Signal length
nb_section=floor(N/(window_label-window_label*overlap_label)); % Number of sections
n_section=1:nb_section; % Section number (ID)
start_time=n_section*(window_label-window_label*overlap_label)-1; % Start time of each section
start_sample=start_time*fn+1; % Start sample of each section
duration_sample=window_label*fn; % Duration of the section in sample

%% POWER BAND OF EACH SECTION
xss_section=reshape(xss, [duration_sample, length(xss)/duration_sample]); % Each section in a column
power_band=bandpower(xss_section, fn, band); % Power band of each section

%% CS and NCS
CS=power_band>threshold; 
NCS=1-CS; 

%% CS REMOVING
xss_section_NCS=xss_section.*NCS; % Each CS column is filled with 0
xss_NCS=reshape(xss_section_NCS, [1 size(xss_section_NCS,1)*size(xss_section_NCS,2)]); % Reshape in a single line
xsc=xss_NCS(xss_NCS~=0); % Signal without CS


end

