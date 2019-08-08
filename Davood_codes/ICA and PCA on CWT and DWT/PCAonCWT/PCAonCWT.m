clear
close all
clc


%%% loading a sample signal
[data,fs]=audioread('NA004.mp3');   
data=downsample(data,8); fs=fs/8;

t1=20;
t2=30;
data=data(t1*fs:t2*fs);
[CW,f_cw]=cwt(data,fs);
[a,b]=find(f_cw>10 & f_cw<1000);
CW=real(CW(a,:));
f_cw=f_cw(a);
j=(1:floor((size(CW,1)/2)))*2; %% size reduction
CW=real(CW(j,:));
f_cw=f_cw(j);
[Xh,U,Y,gamma]=mypca(CW,size(CW,1));



for i=1:size(CW,1)
    CWi=real(CW(i,:))./(3.*std(CW(i,:)));
    subname=['CWT_' int2str(i)];
    filename = ['NA004_' subname '.wav'];
    audiowrite(filename,CWi,fs)
end


for i=1:size(Y,1)
    Yi=real(Y(i,:))./(3.*std(Y(i,:)));
    subname=['PCAonCWT_' int2str(i)];
    filename = ['NA004_' subname '.wav'];
    audiowrite(filename,Yi,fs)
end




figure;
subplot 511
plot((1:size(data,1))./fs,data,'black'); xlim([0 size(data,1)./fs])
subplot 512
plot((1:size(CW,2))./fs,CW(30,:),'r'); xlim([0 size(data,1)./fs])
subplot 513
plot((1:size(Y,2))./fs,Y(30,:),'b'); xlim([0 size(data,1)./fs])
subplot 514
plot((1:size(CW,2))./fs,CW(5,:),'r'); xlim([0 size(data,1)./fs])
subplot 515
plot((1:size(Y,2))./fs,Y(5,:),'b'); xlim([0 size(data,1)./fs])
saveas(gcf,'PcaOnCwt.jpg')

