clear
close all
clc


%%% loading a sample signal
[data,fs]=audioread('NA004.mp3');   
data=downsample(data,8); fs=fs/8;

t1=20;
t2=30;
data=data(t1*fs:t2*fs);
N=6; % no. of levels
data_e1=data(1:end-mod(size(data,1),2^N));  % correcting the length of signal
[SWA,SWD]=swt(data_e1,N,'dmey'); % SWT

%%% normalising the SWT results, by 3*variance
SWA=SWA./(3.*repmat(std(SWA,0,2),1,size(SWA,2)));
SWD=SWD./(3.*repmat(std(SWD,0,2),1,size(SWD,2)));


%%% saving the decomposed data 
for i=1:N
    subname=['SWT level ' int2str(i) ' approximate'];
    filename = ['NA004_' subname '.wav'];
    audiowrite(filename,SWA(i,:),fs)
    subname=['SWT level ' int2str(i) ' detail'];
    filename = ['NA004_' subname '.wav'];
    audiowrite(filename,SWD(i,:),fs)
end





SWC= real([SWA ; SWD]);
B=jadeR(SWC,size(SWC,1));
Y=B*SWC;



for i=1:size(Y,1)
    Yi=Y(i,:)./(3.*std(Y(i,:))); %% uncomment to remove outliers
%     Yi=real(Y(i,:));   %% uncomment to keep outliers
    subname=['ICAonSWT' int2str(i)];
    %title(subname);
    filename = ['NA004_' subname '.wav'];
    audiowrite(filename,Yi,fs)
end

