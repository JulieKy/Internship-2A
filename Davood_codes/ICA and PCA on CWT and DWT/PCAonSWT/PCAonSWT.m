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
SWC= real([SWA ; SWD]);
[Xh,U,Y,gamma]=mypca(SWC,size(SWC,1));




for i=1:size(Y,1)
    Yi=Y(i,:)./(3.*std(Y(i,:))); %% uncomment to remove outliers
%     Yi=real(Y(i,:));   %% uncomment to keep outliers
    subname=['PCAonSWT' int2str(i)];
    filename = ['NA004_' subname '.wav'];
    audiowrite(filename,Yi,fs)
end



for i=1:size(SWC,1)
    Si=SWC(i,:)./(3.*std(SWC(i,:))); %% uncomment to remove outliers
%     Yi=real(Y(i,:));   %% uncomment to keep outliers
    subname=['SWT' int2str(i)];
    filename = ['NA004_' subname '.wav'];
    audiowrite(filename,Si,fs)
end

