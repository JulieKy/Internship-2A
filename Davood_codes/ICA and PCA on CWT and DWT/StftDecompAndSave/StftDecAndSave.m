clear
close all


%%% loading a sample signal
[data,fs]=audioread('nasv\NA004.mp3');   

%%% ploting the sample signal
figure; 
plot((1:size(data,1))/fs, data);   % ploting the sample signal
xlabel('time(sec)'); 
ylabel('scaled amplitude');
axis([0 inf -2 2]);


%%% Calculating and plotting STFT
% spectrogram(data,128,127,[],fs,'yaxis');
% savefig('STFT.fig');

[ST,f_st,t_st]=spectrogram(data,128,127,[],fs,'yaxis');
% figure;
% surface((1:size(ST,2))./fs,f_st,abs(ST));

%%% Saving the results
% save 'tempSTFT.mat'


%%% saving the decomposed signals
N=1;
M=128;
for i=N:M
    STi=real(ST(i,:))./(3.*std(ST(i,:)));
    subname=['STFT freq ' int2str(f_st(i))];
    %title(subname);
    filename = ['NA004_' subname '.wav'];
    audiowrite(filename,STi,fs)
end

