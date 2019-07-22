%% TEST AUTRES FEATURES

clear all; close all;

% Lecture 
[x, Fs]=audioread('1.mp3');

% Preprocessing 
xs=resample(x,4000,Fs);
fn=4000;
y = filterbp(xs,fn);

% Autre
N = length(xs);
time_axis = (1:N)/fn;

%% Linear Predictive Coding (LPC)
[lpc_coeff, var_error] = lpc(y,3);
est_x = filter([0 -lpc_coeff(2:end)],1,xs);

figure,
plot(time_axis,xs); hold on
plot(time_axis,est_x,'--'); hold off
xlabel('Time (s)'); ylabel('Amplitude'); title('LPC?'); legend('Original signal','LPC estimate');

%% Line Spectral Frequencies (LSF)
lsf_coeff = poly2lsf(lpc_coeff); 

%% Brightness : spectrum centroid
% Calculate the centroid of the power spectrum over time. Calculate the centroid for 50 ms Hamming windows of data with 25 ms overlap. Use the range from 62.5 Hz to fs/2 for the centroid calculation. Plot the results.
centroid = spectralCentroid(y,fn, ...
                            'Window',hamming(round(0.05*fn)), ...
                            'OverlapLength',round(0.025*fn), ...
                            'Range',[100,1000]); % verifier le range

t = linspace(0,size(y,1)/fn,size(centroid,1));
figure,
plot(t,centroid);
xlabel('Time (s)'); ylabel('Centroid (Hz)'); title('Spectral centroïd');

% Trouver les picks, les écarts moyens de temps qui les séparent
 [pks,locs] = findpeaks(centroid);
 t_pks=t(locs); 
 mean_tPks=mean(abs(diff(t_pks)));
 
 % Quoi d'autre? La hauteur? 

%% Tonality

%% Spectrogram
[s,f,t] = spectrogram(y); % Quoi faire avec le spectrogramme? Récuperer les zones avec le plus d'intensités en foncyion du temps et de la fréquence? 