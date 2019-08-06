function [threshold] = crying_learning(names_cell)
%CRYING_LEARNING:  Label the crying section, and learn where are the CS by using the Power Ratio Tool

%% INPUTS AND OUTPUTS
%  -- Inputs --
% names_cell:
% -- Outputs --


%% INITIALISATION
% -- Data
observators=3;
samples=37;
end_sample=60; % End of the signal (hypotesis: length of the signal=60s)

% -- Parameters for labelling
window=3;
overlap=25/100;

% -- Parameters for the power ratio
pass_band=[0:2000];
band_width=100;


path = pwd;

%% LABELLING
[label_final, coef_KAPPA]=labelling(observators,samples, end_sample, window, overlap);


%% POWER RATIO

pxx_NCS=[];
band_NCS=[];
PR_NCS=[];

pxx_CS=[];
band_CS=[];
PR_CS=[];

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
    [pxx_NCS_signal, band_NCS_signal, PR_NCS_signal, freq]=power_ratio_band(xss, i, fn, window, overlap, label_final, pass_band, band_width, flag_section);
    pxx_NCS=[pxx_NCS_signal; pxx_NCS];
    band_NCS=[band_NCS_signal; band_NCS];
    PR_NCS=[PR_NCS_signal; PR_NCS];
    f=freq;
    
    % -- CS power ratio
    flag_section=1; % 1 for CS
    [pxx_CS_signal, band_CS_signal, PR_CS_signal, freq]=power_ratio_band(xss, i, fn, window, overlap, label_final, pass_band, band_width, flag_section);
    if (pxx_CS_signal~=0) % CS present in the signal
        pxx_CS=[pxx_CS_signal; pxx_CS];
        band_CS=[band_CS_signal; band_CS];
        PR_CS=[PR_CS_signal; PR_CS];
    end
    
end

% NCS: average on all signals
pxx_NCS_mean=mean(pxx_NCS);
band_NCS_mean=mean(band_NCS);
PR_NCS_mean=mean(PR_NCS);

% CS: average on all signals
pxx_CS_mean=mean(pxx_CS);
band_CS_mean=mean(band_CS);
PR_CS_mean=mean(PR_CS);



%% DISPLAY

% % -- NCS periodogram with means
% figure,
% hax=axes;
% 
% % Display periodogram
% plot(f, pxx_NCS_mean,'LineWidth',2);
% hold on
% 
% band_end=length(f);
% f_interval=length(f)*band_width/(pass_band(end)-pass_band(1));
% 
% % Display lines
% for n_band = 1 : length(f)/f_interval
%     band_start=band_end-(floor(f_interval));
%     line([f(band_start) f(band_start)],get(hax,'YLim'), 'Color',[0 0 0]); % Vertical lines differentiating the frequency bands
%     line([f(band_start) f(band_end)], [band_NCS_mean(n_band) band_NCS_mean(n_band)], 'LineWidth',2, 'Color',[1 0 0]); % Horizontal lines representing the periodogram mean of each frequency band
%     band_end=band_start;
% end
% 
% hold off
% legend('Periodogram', 'Frequency bands', 'Mean')
% title('Welch Periodogram mean for NCS')
% 
% 
% % -- CS periodogram with means
% figure,
% hax=axes;
% 
% band_end=length(f);
% f_interval=length(f)*band_width/(pass_band(end)-pass_band(1));
% 
% 
% % Display periodogram
% plot(f, pxx_CS_mean,'LineWidth',2);
% hold on
% 
% % Display lines
% for n_band = 1 : length(f)/f_interval
%     band_start=band_end-(floor(f_interval));
%     line([f(band_start) f(band_start)],get(hax,'YLim'), 'Color',[0 0 0]); % Vertical lines differentiating the frequency bands
%     line([f(band_start) f(band_end)], [band_CS_mean(n_band) band_CS_mean(n_band)], 'LineWidth',2, 'Color',[1 0 0]); % Horizontal lines representing the periodogram mean of each frequency band
%     band_end=band_start;
% end
% 
% hold off
% legend('Periodogram', 'Frequency bands', 'Mean')
% title('Welch Periodogram mean for CS')

% Median with inter-quartile range
figure(),
plot(f(f~=0) , median(pxx_NCS(:,(f~=0))), 'Color', [1 0 0]); hold on
plot(f(f~=0) , median(pxx_CS(:,(f~=0))), 'Color', [0 0 1]);
plot(f(f~=0) , prctile(pxx_NCS(:,(f~=0)),25), 'LineStyle', '--', 'Color', [1 0.4 0.4]); 
plot(f(f~=0) , prctile(pxx_CS(:,(f~=0)),25), 'LineStyle', '--', 'Color', [0.4 0.4 1]); 
plot(f(f~=0) , prctile(pxx_NCS(:,(f~=0)),75), 'LineStyle', '--', 'Color', [1 0.4 0.4]); 
plot(f(f~=0) , prctile(pxx_CS(:,(f~=0)),75), 'LineStyle', '--', 'Color', [0.4 0.4 1]); 

xlabel('Frequency (in Hz)');
ylabel('Power Spectrum');
title('Average Power Spectrum for NCS and CS');
legend('Average power spectrum NCS',  'Average power spectrum CS', 'Interquartile range NCS','Interquartile range CS');
hold off


threshold=1;
end
