clear
close all
clc


%%% loading a sample signal
[data,fs]=audioread('NA004.mp3');   
data=downsample(data,8); fs=fs/8;

t1=20;
t2=30;
M=2;N=129;
data=data(t1*fs:t2*fs)';
audiowrite('NA004forProcess.wav',data,fs)

%%% Calculating STFT
[ST,f_st,t_st]=spectrogram(data,128,127,[],fs,'yaxis');
ST=real(ST);
%%% PCA ...
STforJADE=ST(M:N,:);
j=(1:floor((size(STforJADE,1))/4))*4; %% size reduction
STforJADE=STforJADE(j,:);
B=jadeR(STforJADE,size(STforJADE,1));
Y=B*STforJADE;
Y=real(Y);

% 
% for i=1:size(STforJADE,1)
%     STi=STforJADE(i,:)./(3.*std(STforJADE(i,:))); %% uncomment to remove outliers
% %     STi=real(ST(i,:));   %% uncomment to keep outliers
%     subname=['STFT_' int2str(i)];
%     filename = ['NA004_' subname '.wav'];
%     audiowrite(filename,STi,fs)
% end
% 
% 
% 
% for i=1:size(Y,1)
%     Yi=Y(i,:)./(3.*std(Y(i,:))); %% uncomment to remove outliers
% %     Yi=real(Y(i,:));   %% uncomment to keep outlier
%     subname=['JADEonSTFT' int2str(i)];
%     filename = ['NA004_' subname '.wav'];
%     audiowrite(filename,Yi,fs)
% end



figure;
subplot 511
plot((1:size(data,2))./fs,data,'black'); xlim([0 size(data,2)./fs]);
subplot 512
plot((1:size(STforJADE,2))./fs,STforJADE(21,:),'r');
subplot 513
plot((1:size(Y,2))./fs,Y(18,:),'b');
saveas(gcf,'JadeOnStft.jpg')
subplot 514
plot((1:size(STforJADE,2))./fs,STforJADE(5,:),'r');
subplot 515
plot((1:size(Y,2))./fs,Y(8,:),'b');
saveas(gcf,'JadeOnStft.jpg')

audiowrite('JadeBestLung.wav',Y(18,:)./(3.*std(Y(18,:))),fs)
audiowrite('StftBestLung.wav',STforJADE(21,:)./(3.*std(STforJADE(21,:))),fs)

audiowrite('JadeBesHeart.wav',Y(8,:)./(3.*std(Y(8,:))),fs)
audiowrite('StftBestHeart.wav',STforJADE(5,:)./(3.*std(STforJADE(5,:))),fs)



