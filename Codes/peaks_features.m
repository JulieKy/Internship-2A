function [nb_pks,  f_higherPk, dif_higherPks] = peaks_features(pxx, pxx_smooth, f)
%% INPUT AND OUTPUT
% -- Inputs --
% pxx: periodogram 
% pxx_smooth: smoothed or fitted
% f: frequency
% -- Outputs --
% nb_pks: number of peaks
% f_pks: frequency of these peaks
% dif_higherPks: 2 higher peaks frequency differences

%% NUMBER OF PEAKS
 [pks,locs] = findpeaks(pxx_smooth);
 nb_pks=length(pks);

 %% FREQUENCY OF PEAKS (in ascending order im term of PSD)
 [pksOrder,order] = sort(pks); % Rank the peaks in ascending order
 f_pks=f(locs(order)); 
 
 %%  2 HIGHER PEAKS FREQUENCY 
 nb_higherPks=min(length(pks),2);
 higherPks=[]; 
 argmax_higherPks=[]; 
 f_higherPks=[];
 for i=0:nb_higherPks-1
     higherPks=[higherPks, pksOrder(end-i)];
     argmax_higherPks=[argmax_higherPks, order(end-i)];
     f_higherPks=[f_higherPks, f(locs(argmax_higherPks(i+1)))];
 end
 f_higherPk=f_higherPks(1);
 dif_higherPks=f_higherPks(1)-f_higherPks(end); % 0 if there is only one peak
 
 
%  %%  DISPLAY 
% figure,plot(f,pxx); hold on
% plot(f,pxx_smooth,f(locs),pks,'or','MarkerEdgeColor','green'); 
% plot(f_higherPks,higherPks,'r*','MarkerEdgeColor','red'); hold off
% xlabel('Frequency (0-1000Hz)'),ylabel('PSD (dB/Hz)'), title('Periodogram Welch','fontsize',14,'interpreter','latex'), legend('pxx', 'pxx\_smooth','interpreter','latex');

end

