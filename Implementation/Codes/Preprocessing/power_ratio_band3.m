function [PR]=power_ratio_band3( pass_band, fn, pure_sections, band_width)
%POWER RATIO:  Calculate the power ratio of NCS/CS (depending on flag_section) of a signal

%% INPUTS AND OUTPUTS
%  -- Inputs --
% xss: Signal without treatment on CS/NCS
% signal_n: Number of the signal
% fn: Sampling frequency
% window: Window used for labelling
% overlap: Overlap used for labelling
% label_final: Annotated labels of the signal bank
% flag_section: 1=CS; 0=NCS
% -- Outputs --

%% INITIALISATION
sections_nb=size(pure_sections, 2); % Number of sections
band_nb=(pass_band(end)-pass_band(1))/band_width;
PR_sections=zeros(band_nb, sections_nb); % Power Ratio for different frequency bands, on all pure sections

%% POWER RATIO
band_width=1;
% -- For each band, computation of the PR 
for b=pass_band(1)+1:band_width:pass_band(end)+1
    band=[b, b+band_width]; % Frequency band
    band=[0 2000];
    PR_sections(:,b)=bandpower(pure_sections, fn, band); % Computation of the PR for each section
end

%% OUTPUT
PR=mean(PR_sections); % Mean of the PR of each section, on the different bands
    
end


