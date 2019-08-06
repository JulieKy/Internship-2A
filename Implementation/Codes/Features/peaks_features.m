function [nb_pks,  f_higherPk, dif_higherPks] = peaks_features(y, f, nb_higherPks_method, tag, t)
%% INPUT AND OUTPUT
% -- Inputs --
% y: periodogram or spectrogram
% f: frequency
% nb_higherPks_method: number of higher peaks wanted 
% tag: Type of periodogram or spectrogram
% t: time (useful for the spectrogram)
% -- Outputs --
% nb_pks: number of peaks
% f_pks: frequency of these peaks
% dif_higherPks: 2 higher peaks frequency differences


%% NUMBER OF PEAKS
 [pks,locs] = findpeaks(y);
 nb_pks=length(pks);

 %% FREQUENCY OF PEAKS (in ascending order im term of PSD)
 [pksOrder,order] = sort(pks); % Rank the peaks in ascending order
 f_pks=f(locs(order));    
 
 %%  HIGHER PEAKS FREQUENCY 
 nb_higherPks=min(length(pks),nb_higherPks_method);
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
 
 
%%  HIGHER PEAKS TIME
 
 %%  DISPLAY 

if strcmp('spectrogram', tag)==0
    
    % Display the periodogram
    hold on
    plot(f,y,f(locs),pks,'or','MarkerEdgeColor','green');  
    plot(f_higherPks,higherPks,'r*','MarkerEdgeColor','red'); hold off
    xlabel('Frequency (0-1000Hz)'),ylabel('PSD (dB/Hz)');

    if strcmp('periodogram_MAF', tag)==1
        title('Welch Periodogram with Moving Average Filter (MAF)','fontsize',14,'interpreter','latex'),
        legend('pxx', 'pxx\_smooth with MAF','interpreter','latex');
    elseif strcmp('periodogram_GMM', tag)==1
        title('Welch Periodogram with Gaussian Mixture Model (GMM)','fontsize',14,'interpreter','latex'),
        legend('pxx', 'pxx\_smooth with GMM','interpreter','latex');
    end
    
else 
     % Display the spectrogram
    print('prouuttt')
end 


end
