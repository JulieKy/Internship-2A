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
j=(1:floor((size(CW,1)-44)/2))*2; %% size reduction
CWforJADE=real(CW(j,:));

B=jadeR(CWforJADE,size(CWforJADE,1));
Y=B*CWforJADE;




for i=1:size(Y,1)
    Yi=real(Y(i,:))./(5.*std(Y(i,:)));
    subname=['JADEonCWT_' int2str(i)];
    filename = ['NA004_' subname '.wav'];
    audiowrite(filename,Yi,fs)
end



for i=1:size(CWforJADE,1)
    CWforJADEi=real(CWforJADE(i,:))./(5.*std(CWforJADE(i,:)));
    subname=['CWTforJADE_' int2str(i)];
    filename = ['NA004_' subname '.wav'];
    audiowrite(filename,CWforJADEi,fs)
end




