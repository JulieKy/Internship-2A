% This function computes spectral features 
%
function [output_spectral_features,pxx,f,foct,spower,I,S] = spectral_features(x,fn) % See Fae's comments
%% INPUT AND OUTPUT
% input = x 'audio signal'
% Fs 'sampling frequency
% outputs =  spectral features including
% meanPSD: mean frequency of power spectrum
% stdPSD: std frequency of power spectrum
% medPSD:  median frequency of power spectrum
% bw: 3db bandwidth
% p25: 1st quartile of power spectrum
% p75: 3rd quartile of power spectrum
% IQR: inter quartile range
% TP: total power in 100-1000 Hz
% p100_200: power ratio: 100-200 hz/TP
% p200_400: power ratio: 200-400 hz/TP
% p400_800: power ratio: 400-800 hz/TP
% spectrum_slope2: spectrum slope
% r_square2: R^2 statistics (linear regression for slope)
% nbPeaks: Numbers of peaks in the Welch Periodogram (Julie)
% DifFreqFitPeaks: Frequency difference between the two higher peaks (Julie)
% DifHighFitPeaks: Amplitude difference between the two higher peaks (Julie)


%% PERIODOGRAM WELCH (MATHIEU)

% -- Implementation
wind = ones(1,floor(0.125*fn)); % 125 ms
nover = floor(length(wind)/4); % 25% overlap
nfft = 2^(nextpow2(length(wind))-1); % nfft
[pxx,f] = pwelch(x,wind,nover,nfft,fn); % Welch
f = f(f<=1000); % [0-800Hz] band of interest % Fae changed it to 1000Hz
pxx=pxx(1:length(f));
power = 1                  
xdata2 = f; % used for linear regression
ydata2 = power; % used for linear regression


%% PERIODOGRAM WITH SMOOTH and FIT METHODS (Julie)

%% Moving Average Filter
% -- Numbers of peaks
pxx_smooth = smooth(pxx);
[pks,locs] = findpeaks(pxx_smooth);
nbPeaks_MAF=length(pks);

% % -- Display figure
% figure,plot(f,pxx); hold on
% plot(f,pxx_smooth,f(locs),pks,'or'); hold off
% xlabel('Frequency (0-1000Hz)','fontsize',11,'interpreter','latex'),ylabel('PSD (dB/Hz)','fontsize',11,'interpreter','latex'), title('Periodogram Welch','fontsize',14,'interpreter','latex'), legend('pxx', 'pxx\_smooth','interpreter','latex');


%% Gaussian Mixture Model
% -- Fitting Gaussians to the spectrum
 fi = fit(f,pxx,'gauss4');
 fi1=fi.a1*exp(-((f-fi.b1)./fi.c1).^2);
 fi2=fi.a2*exp(-((f-fi.b2)./fi.c2).^2);
 fi3=fi.a3*exp(-((f-fi.b3)./fi.c3).^2);
 fi4=fi.a4*exp(-((f-fi.b4)./fi.c4).^2);
 fi_tot=fi1+fi2+fi3+fi4;
 
 % -- Gaussian parameters
a=[fi.a1, fi.a2, fi.a3, fi.a4];
b=[fi.b1, fi.b2, fi.b3, fi.b4];
c=[fi.c1, fi.c2, fi.c3, fi.c4];
 
% -- Number of peaks
 [pks_GMM,locs_GMM] = findpeaks(fi_tot);
 nbPeaks_GMM=length(pks_GMM);
 
 % -- Frequency of peaks (in ascending order im term of PSD)
 [pksOrder_GMM,order_GMM] = sort(pks_GMM); % Rank the peaks in ascending order
 f_pks_GMM=f(locs_GMM(order_GMM)); 
 
 % -- 2 higher peaks frequency differences
 nb_higher_pks=min(length(pks_GMM),2);
 higher_pks=[]; 
 argmax_higher_pks=[]; 
 f_higher_pks=[];
 for i=0:nb_higher_pks-1
     higher_pks=[higher_pks, pksOrder_GMM(end-i)];
     argmax_higher_pks=[f_higher_pks, order_GMM(end-i)];
     f_higher_pks=[f_higher_pks, f(locs_GMM(argmax_higher_pks(i+1)))];
 end
 dif=f_higher_pks(1)-f_higher_pks(end); 
  
% % -- Display figure
figure,
plot(f,pxx); hold on
plot(f,fi1); plot(f,fi2); plot(f,fi3); plot(f,fi4); plot(f,fi_tot); % Gaussian fit
plot(f_higher_pks,higher_pks,'or'); hold off % 2 higher peaks
xlabel('Frequency (0-1000Hz)','fontsize',11,'interpreter','latex'),ylabel('PSD (dB/Hz)','fontsize',11,'interpreter','latex'), title('Periodogram Welch and Gaussians fitting','fontsize',14,'interpreter','latex'), legend('periodogram','fi1','fi2','fi3','fi4','fi_tot','interpreter','latex');



%% PERIODOGRAM ANALYSIS (Notes by Fae)
% Note that in the Peruvian paper it mentions: "the slope of the linear regression line, fit to
% spectrum P in logarithmic axes. The power spectrum, when plotted in dB as
% 20*log (P/Pmin)". Let's assume that Pmin is the power at 1000Hz, then the
% 'relative' power in db is: (however I tink there was a mistake in the article and it should be 10*log10(p/pmin))
spower=10*log10(pxx/pxx(end));

%% spectrum parameters
% some statistical spectral parameters:
meanPSD=meanfreq(pxx,f); %  The mean frequency of a power spectrum
stdPSD=sqrt(sum(pxx.*((f-meanPSD).^2))/sum(pxx));%  The std of a power spectrum
medPSD=medfreq(pxx,f); %  The median frequency of a power spectrum
bw = powerbw(pxx,f); % 3db bandwidth of power spectrum
p25=percentilefreq(pxx,f,25); % 25 percentile (1st quartile) of spectral power
p75=percentilefreq(pxx,f,75); % 75 percentile (3rd quartile) of spectral power
IQR=p75-p25; % Interquartile range of spectrum
 % power ratios of various frequency (Hz) bands 
 TP=bandpower(pxx,f,[100 f(end)],'psd');% total power for 100-1000Hz
p100_200 = bandpower(pxx,f,[100 200],'psd')/TP;
p200_400 = bandpower(pxx,f,[200 400],'psd')/TP;
p400_800 = bandpower(pxx,f,[400 800],'psd')/TP;

%% Slope of the regression line (SL) in db/octave (Notes by Fae)
% for octave we need should convert the frequency like this (taking the base as mean frequency):  
foct=log2(f/meanPSD);
% To fit a linear reg line:
mdl = fitlm(foct(foct>0),spower(foct>0));
I=mdl.Coefficients{1,1};% Intercept
S=mdl.Coefficients{2,1};% slope
limOct=log2(1000/meanPSD);
% figure,
% plot(foct,spower,'k*:'),xlim([0,limOct]),title('Spectrum Slope'),xlabel('Hz with oct divisions'),ylabel('dB')
% hold on
% plot(foct,I+S*foct)
% hold off 
% set(gca,'XTick',[0:ceil(limOct)] );
% set(gca,'XTickLabel', [0 2.^[1:ceil(limOct)]].*round(meanPSD) );


spectrum_slope2 =S;
r_square2=mdl.Rsquared.Adjusted;


%% RESULT

%combined final output

output_spectral_features = [meanPSD;stdPSD;medPSD;bw;p25;p75;IQR;TP;p100_200;p200_400;p400_800;spectrum_slope2;r_square2;nbPeaks;DifFreqFitPeaks;DifHighFitPeaks];

end

