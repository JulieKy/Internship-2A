clear
close all

%%% loading a sample signal
[data,fs]=audioread('nasv\NA004.mp3');

%%% ploting the sample signal
figure;
plot((1:size(data,1))/fs, data);
xlabel('time(sec)'); 
ylabel('scaled amplitude');
axis([0 inf -2 2]);

%%% Calculating and plotting CWT
[CW,f_cw]=cwt(data,fs);
% figure;
% surface((1:size(CW,2))./fs,f_cw,abs(CW),,'EdgeColor','none');

%%% saving the results
% save('temp.mat','CW','-v7.3')
% save('temp.mat','f_cw','-append');




%%% saving the decomposed signals
N=1;
M=170;
for i=N:M
    CWi=real(CW(i,:))./(3.*std(CW(i,:)));
    subname=['CWT freq ' int2str(f_cw(i))];
    filename = ['NA004_' subname '.wav'];
    audiowrite(filename,CWi,fs)
end
    