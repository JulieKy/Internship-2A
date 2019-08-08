% This function computes spectral features 
%
function [output_spectral_features,pxx,f,foct,spower,I,S] = spectral_features_4p(x,Fs) % See Fae's comments
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
% 


%% BAND PASSFILTER (FATEMA)
%filtering the input signal
%input should be 'minimum_seg_longest_noncrying.wav' created in
%noncrying_samples.m function
%since spectral features are required of the filtered non crying segmnets
%if one wishes to conduct on raw signal or any other, the first 4 lines
%should be commented and fn = Fs should be written


%%% the 4 nex lines are commented by davood
% filterbp(x,Fs);
% [filtered_signal,fsfilteredsignal]=audioread('filtered_signal.wav');
% x = filtered_signal;
% fn = fsfilteredsignal;
fn=Fs;

%% PERIODOGRAM WELCH (MATHIEU)

wind = ones(1,floor(0.125*fn)); % 125 ms
nover = floor(length(wind)/4); % 25% overlap
nfft = 2^(nextpow2(length(wind))-1); % nfft
[pxx,f] = pwelch(x,wind,nover,nfft,fn); % Welch
f = f(f<=625); % [0-800Hz] band of interest % Davood changed it to 625 Hz
pxx=pxx(1:length(f));
power = 10*log10(pxx); % computation of power in dB for band of interest

xdata2 = f; % used for linear regression
ydata2 = power; % used for linear regression

% % Display figure
% figure,plot(f,pxx);
% xlabel('Frequency (0-800Hz)'),ylabel('PSD (dB/Hz)'), title('Periodogram Welch');



%% Gaussian fit
% fitting Gaussians to the spectrum, not using it for now, maybe later
 fi = fit(f,pxx,'gauss2');
 fi1=fi.a1*exp(-((f-fi.b1)./fi.c1).^2);
 fi2=fi.a2*exp(-((f-fi.b2)./fi.c2).^2);
 fi3=fi1+fi2;
% figure,
% plot(f,pxx);hold on
% plot(f,fi1)
% plot(f,fi2)
% plot(f,fi3)


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
 TP=bandpower(pxx,f,[125 625],'psd');% total power for 100-1000Hz
p125_250 = bandpower(pxx,f,[125 250],'psd')/TP;
p250_375 = bandpower(pxx,f,[250 375],'psd')/TP;
p375_500 = bandpower(pxx,f,[375 500],'psd')/TP;
p500_625 = bandpower(pxx,f,[500 625],'psd')/TP;

%% Slope of the regression line (SL) in db/octave (Notes by Fae)
% for octave we need should convert the frequency like this (taking the base as mean frequemcy):  
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
% set(gca,'XTick',[0:ceil(limOct)] );
% set(gca,'XTickLabel', [0 2.^[1:ceil(limOct)]].*round(meanPSD) );

spectrum_slope2 =S;
r_square2=mdl.Rsquared.Adjusted;


%% RESULT

%combined final output

output_spectral_features = [meanPSD;stdPSD;medPSD;bw;p25;p75;IQR;TP;p125_250;p250_375;p375_500;p500_625;spectrum_slope2;r_square2];

end

