%%%%%%%%%%%%%%%%%%%
%% Initializing
clear
close all
clc

files=dir('Samples\*.mp3'); %% Defining address of the folder containing the audio recorde



fcl=150;  %% Low cutoff frequency
fch=800;    %% High cutoff frequency
OutputSheet=cell(size(files,1)+1,7);
OutputSheet(1,1:7)={'name','LPCC1', 'LPCC2', 'LPCC3', 'LPCC4', 'LPCC5', 'LPCC6'};
%% Doing for each of the records ... 
for i=1:size(files,1)
    [data,fs]=audioread(['Samples\' files(i).name]);    %% read the record
   
    %%
    % Bandpass filtering 1 :     
    dlp=fdesign.lowpass('Fp,Fst,Ap,Ast', .9*fch, 1.1*fch, .1, 60,fs);
    Hdlp = design(dlp,'equiripple');
    % fvtool(Hdlp)
    % measure(Hd)

    dhp = fdesign.highpass('Fst,Fp,Ast,Ap', .9*fcl, 1.1*fcl, 60, .1, fs);
    Hdhp = design(dhp,'equiripple');
    % fvtool(Hdhp)
    % measure(Hdhp)

    data=filter(Hdlp,(filter(Hdhp,data)));

    %% or Bandpass filtering 2 :
%     d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',.9*fcl,1.1*fcl,0.9*fch,1.1*fch,60,5,60,fs);
%     Hd = design(d,'equiripple');
%     data=filter(Hd,data);
    %% LPCC Feature extracting
    [LPCC,g] = lpc(data,6);
    OutputSheet(i+1,1)=cellstr(files(i).name);
    OutputSheet(i+1,2:7)=num2cell(LPCC(2:end));
end

writecell(OutputSheet,'OutputSheet_LPCC_150_800.xls');

