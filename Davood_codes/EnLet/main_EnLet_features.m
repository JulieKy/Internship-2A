%%%%%%%%%%%%%%%%%%%
%% Initializing
clear
close all
clc

files=dir('Samples\*.mp3'); %% Defining address of the folder containing the audio recorde



fcl=150;  %% Low cutoff frequency
fch=800;    %% High cutoff frequency
OutputSheet=cell(size(files,1)+1,7);
OutputSheet(1,1:13)={'name','EnApp1', 'EnDet1', 'EnApp2', 'EnDet2', 'EnApp3', 'EnDet3', 'EnApp4' 'EnDet4', 'EnApp5', 'EnDet5', 'EnApp6', 'EnDet6'};
%% Doing for each of the records ... 
for i=1:size(files,1)
    [data,fs]=audioread(['Samples\' files(i).name]);    %% read the record
    data=downsample(data,8);
    fs=fs/8;
    
    approx=data;
%     dim=floor(fs/100);
    dim=fs/500;
    N=6;
    for j=1:N
       disp([num2str(j) ' of ' num2str(N) ' levels, ' num2str(i) ' of ' num2str(size(files,1)) ' cases, ...']);
       [approx, det]=dwt(approx,'dmey');
       approxapen=ApEn(dim,std(approx),approx,1);
       detapen=ApEn(dim,std(det),det,1);
       ApenOnDWT(:,j)=[approxapen;detapen];
       dim=floor(dim/2);
    end
    OutputSheet(i+1,1)=cellstr(files(i).name);
    OutputSheet(i+1,2:13)=num2cell(ApenOnDWT(:)');
end

writecell(OutputSheet,'OutputSheet_EnLet.xls');

