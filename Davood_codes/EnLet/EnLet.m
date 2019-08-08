clear
close all
clc
N=6;

%%% loading a sample signal
[data,fs]=audioread('NA004.mp3'); 
approx=data;
dim=fs/10;
for i=1:N
    [approx, det]=dwt(approx,'dmey');
    approxapen=ApEn(dim,.2*std(approx),approx,1);
    detapen=ApEn(dim,.2*std(det),det,1);
    ApenOnDWT(:,i)=[approxapen;detapen];
    dim=dim/2;
end
    
    