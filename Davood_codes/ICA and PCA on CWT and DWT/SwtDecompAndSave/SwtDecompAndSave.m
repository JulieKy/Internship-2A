clear
close all

%%% loading a sample signal
[data,fs]=audioread('nasv\NA004.mp3');

%%% ploting the sample signal
figure;
%subplot(411)
plot((1:size(data,1))/fs, data);
xlabel('time(sec)'); 
ylabel('scaled amplitude');
axis([0 inf -2 2]);


%%% Calculating and plotting SWT
N=6; % no. of levels
data_e1=data(1:end-mod(size(data,1),2^N));  % correcting the length of signal
[SWA,SWD]=swt(data_e1,N,'dmey'); % SWT

%%% normalising the SWT results, by 3*variance
SWA=SWA./(3.*repmat(std(SWA,0,2),1,size(SWA,2)));
SWD=SWD./(3.*repmat(std(SWD,0,2),1,size(SWD,2)));

%%% plotting and saving the decomposed data 
figure;
for i=1:N
    subplot(2*N,1,2*i-1); plot((1:size(data_e1,1))/fs,SWA(i,:)); axis([0 inf -2 2]);
    %xlabel('time(sec)'); 
    subname=['SWT level ' int2str(i) ' approximate'];
    %title(subname);
    filename = ['NA004_' subname '.wav'];
    audiowrite(filename,SWA(i,:),fs)
    subplot(2*N,1,2*i); plot((1:size(data_e1,1))/fs,SWD(i,:)); axis([0 inf -2 2]);
    %xlabel('time(sec)');
    subname=['SWT level ' int2str(i) ' detail'];
    %title(subname);
    filename = ['NA004_' subname '.wav'];
    audiowrite(filename,SWD(i,:),fs)
end
