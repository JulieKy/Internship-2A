function [truc] = crying_learning(names_cell)
%CRYING_LEARNING:  Label the crying section, and learn where are the CS by using the Power Ratio Tool

%% INPUTS AND OUTPUTS





%% INITIALISATION
% -- Data
observators=2;
samples=37;
end_sample=60; % End of the signal (hypotesis: length of the signal=60s)

% -- Parameters for labelling
window=3;
overlap=25/100;

% -- Parameters for the power ratio
pass_band=[0:1000];
band_width=100;

%% LABELLING
[label_final, coef_KAPPA]=labelling(observators,samples, end_sample, window, overlap);


%% POWER RATIO

pxx_NCS=[];
band_NCS=[];
PR_NCS=[];

pxx_CS=[];
band_CS=[];
PR_CS=[];


path = pwd;

% For every signal
for i = 1:length(names_cell)
    
    % -- Reading the signals
    tempName=names_cell{i};
    [x,Fs]= audioread([path,'\..\Data\Samples_Belle\',tempName]);
    
    % -- Resampling to 4000 Hz
    xs=resample(x,4000,Fs);
    fn=4000;
    
    % -- Shorten to 60s
    time_sample=60;
    xss=xs(1:time_sample*fn,1);
    
    % -- NCS power ratio
    flag_section=0; % 0 for NCS
    [pxx_NCS_signal, band_NCS_signal, PR_NCS_signal, f]=power_ratio_band(xss, fn, window, overlap, label_final, pass_band, band_width, flag_section);
    pxx_NCS=[pxx_NCS_signal; pxx_NCS];
    band_NCS=[band_NCS_signal; band_NCS];
    PR_NCS=[PR_NCS_signal; PR_NCS];
    
    % -- CS power ratio
    flag_section=1; % 1 for CS
    [pxx_CS_signal, band_CS_signal, PR_CS_signal, f]=power_ratio_band(xss, fn, window, overlap, label_final, pass_band, band_width, flag_section);
    pxx_CS=[pxx_CS_signal; pxx_CS];
    band_CS=[band_CS_signal; band_CS];
    PR_CS=[PR_CS_signal; PR_CS];
    
end

% NCS
pxx_NCS_mean=mean(pxx_NCS);
band_NCS_mean=mean(band_NCS);
PR_NCS_mean=mean(PR_NCS);

% CS
pxx_CS_mean=mean(pxx_CS);
band_CS_mean=mean(band_CS);
PR_CS_mean=mean(PR_CS);


%% DISPLAY

% -- NCS periodogram with means
figure,
hax=axes;

% Display periodogram
plot(f, pxx_NCS_mean,'LineWidth',2);
hold on

band_end=length(f);
f_interval=length(f)*band_width/(pass_band(end)-pass_band(1))

% Display lines
for n_band = 1 : length(f)/f_interval
    band_start=band_end-(floor(f_interval));
    line([f(band_start) f(band_start)],get(hax,'YLim'), 'Color',[0 0 0]); % Vertical lines differentiating the frequency bands
    line([f(band_start) f(band_end)], [band_NCS_mean(n_band) band_NCS_mean(n_band)], 'LineWidth',2, 'Color',[1 0 0]);
    band_end=band_start;
end

hold off
legend('Periodogram', 'Frequency bands', 'Mean')
title('Welch Periodogram mean for NCS')

truc=1;
end

