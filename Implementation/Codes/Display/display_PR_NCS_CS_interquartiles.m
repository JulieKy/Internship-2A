function [] = display_PR_NCS_CS_interquartiles(f,pxx_NCS, pxx_CS, pxx_NCS_mean, pxx_CS_mean, band_width, pass_band, band_NCS_mean, band_CS_mean)
%DISPLAY_PR_NCS_CS_INTERQUARTILES: Display the periodograms of annotated NCS and CS

%% INPUTS AND OUTPUTS
%  -- Inputs --
% f: Vectors of frequencies 
% pxx_NCS: Periodograms of all NCS 
% pxx_CS: Periodograms of all CS 
% pxx_NCS_mean: Mean of NCS periodograms 
% pxx_CS_mean: Mean of CS periodograms
% band_width: Width of the frequency bands
% pass_band: Pass band
% band_NCS_mean: Width of the NCS frequency bands
% band_CS_mean: Width of the CS frequency bands

% %% NCS periodogram with means of each frequency bands
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
% %% CS periodogram with means of each frequency bands
% band_end=length(f);
% 
% figure, 
% hax=axes;
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


%% Median with inter-quartile range
figure,
plot(f(f~=0) , median(pxx_NCS(:,(f~=0))), 'Color', [1 0 0]); hold on
plot(f(f~=0) , median(pxx_CS(:,(f~=0))), 'Color', [0 0 1]);
plot(f(f~=0) , prctile(pxx_NCS(:,(f~=0)),25), 'LineStyle', '--', 'Color', [1 0.4 0.4]); 
plot(f(f~=0) , prctile(pxx_CS(:,(f~=0)),25), 'LineStyle', '--', 'Color', [0.4 0.4 1]); 
plot(f(f~=0) , prctile(pxx_NCS(:,(f~=0)),75), 'LineStyle', '--', 'Color', [1 0.4 0.4]); 
plot(f(f~=0) , prctile(pxx_CS(:,(f~=0)),75), 'LineStyle', '--', 'Color', [0.4 0.4 1]); 

xlabel('Frequency (in Hz)');
ylabel('Power Spectrum');
title('Average Power Spectrum for NCS and CS (labelling: window=1s)', 'FontSize', 14);
legend('Average power spectrum NCS',  'Average power spectrum CS', 'Interquartile range NCS','Interquartile range CS');
hold off

end

