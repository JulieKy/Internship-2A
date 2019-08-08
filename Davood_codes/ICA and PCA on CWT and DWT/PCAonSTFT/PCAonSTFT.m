clear
close all
clc


%%% loading a sample signal
[data,fs]=audioread('NA004.mp3');   
data=downsample(data,4); fs=fs/4;

t1=20;
t2=30;
M=2;N=129;
data=data(t1*fs:t2*fs)';
audiowrite('NA004forProcess.wav',data,fs)

%%% Calculating STFT
[ST,f_st,t_st]=spectrogram(data,128,127,[],fs,'yaxis');
ST=real(ST);
%%% PCA ...
STforPCA=ST(M:N,:)';
[Xh,U,Y,gamma]=mypca(STforPCA',N-M+1);
Y=real(Y);


% 
% for i=M:N
%     STi=ST(i,:)./(3.*std(ST(i,:))); %% uncomment to remove outliers
% %     STi=real(ST(i,:));   %% uncomment to keep outliers
%     subname=['STFT freq ' int2str(f_st(i))];
%     %title(subname);
%     filename = ['NA004_' subname '.wav'];
%     audiowrite(filename,STi,fs)
% end
% 
% 
% 
% for i=1:N-M+1
%     Yi=Y(i,:)./(3.*std(Y(i,:))); %% uncomment to remove outliers
% %     Yi=real(Y(i,:));   %% uncomment to keep outliers
%     subname=['PCAonSTFT' int2str(i)];
%     %title(subname);
%     filename = ['NA004_' subname '.wav'];
%     audiowrite(filename,Yi,fs)
% end



figure;
subplot 511
plot((1:size(data,2))./fs,data,'black'); xlim([0 size(data,2)./fs]);
subplot 512
plot((1:size(STforPCA,1))./fs,STforPCA(:,50),'r');
subplot 513
plot((1:size(Y,2))./fs,Y(50,:),'b');
subplot 514
plot((1:size(STforPCA,1))./fs,STforPCA(:,6),'r');
subplot 515
plot((1:size(Y,2))./fs,Y(6,:),'b');
saveas(gcf,'PcaOnStft.jpg')

audiowrite('PCABestLung.wav',Y(50,:)./(3.*std(Y(50,:))),fs)
audiowrite('StftBestLung.wav',STforPCA(:,50)./(3.*std(STforPCA(:,50))),fs)

audiowrite('PCABesHeart.wav',Y(6,:)./(3.*std(Y(6,:))),fs)
audiowrite('StftBestHeart.wav',STforPCA(:,6)./(3.*std(STforPCA(:,6))),fs)



