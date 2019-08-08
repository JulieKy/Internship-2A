%%%%%%%%%%%%%%%%%%%
%% Initializing
clear
close all
clc

files=dir('Samples\*.mp3'); %% Defining address of the folder containing the audio recorde



fcl=100;  %% Low cutoff frequency
fch=1000;    %% High cutoff frequency
OutputSheet=cell(size(files,1)+1,20);
OutputSheet(1,1:26)={'name','meanPSD','stdPSD','medPSD','bw','p25','p75','IQR','TP','p100_200','p200_400','p400_800','spectrum_slope2','r_square2', 'mfcc1_mean', 'mfcc2_mean', 'mfcc3_mean', 'mfcc4_mean', 'mfcc5_mean', 'mfcc6_mean'};
%% Doing for each of the records ... 
for i=1:size(files,1)
    [data,fs]=audioread(['Samples\' files(i).name]);    %% read the record
    
    %% MFCC 
    mfcc_mean = mean(mfcc(data,fs, 'NumCoeffs', 6,'LogEnergy','Ignore'));
   
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
    %% Feature extracting
    [output_spectral_features,pxx,f,foct,spower,I,S] = spectral_features_3p(data,fs);
    OutputSheet(i+1,1)=cellstr(files(i).name);
    OutputSheet(i+1,2:14)=num2cell(output_spectral_features');
    OutputSheet(i+1,15:20)=num2cell(mfcc_mean);

end

writecell(OutputSheet,'OutputSheet_DavBandpass_100_1000.xls');




OutputSheet=cell(size(files,1)+1,20);
OutputSheet(1,1:20)={'name','meanPSD','stdPSD','medPSD','bw','p25','p75','IQR','TP','p100_200','p200_400','p400_800','spectrum_slope2','r_square2', 'mfcc1_mean', 'mfcc2_mean', 'mfcc3_mean', 'mfcc4_mean', 'mfcc5_mean', 'mfcc6_mean'};
%% Doing for each of the records ... 
for i=1:size(files,1)
    [data,fs]=audioread(['Samples\' files(i).name]);    %% read the record
    
    %% MFCC 
    mfcc_mean = mean(mfcc(data,fs, 'NumCoeffs', 6,'LogEnergy','Ignore'));

%     % Bandpass filtering 1 :     
%     dlp=fdesign.lowpass('Fp,Fst,Ap,Ast', .9*fch, 1.1*fch, 5, 60,fs);
%     Hdlp = design(dlp,'equiripple');
%     % fvtool(Hdlp)
%     % measure(Hd)
% 
%     dhp = fdesign.highpass('Fst,Fp,Ast,Ap', .9*fcl, 1.1*fcl, 60, 5, fs);
%     Hdhp = design(dhp,'equiripple');
%     % fvtool(Hdhp)
%     % measure(Hdhp)
% 
%     data=filter(Hdlp,(filter(Hdhp,data)));

    %% or Bandpass filtering 2 :
%     d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',.9*fcl,1.1*fcl,0.9*fch,1.1*fch,60,5,60,fs);
%     Hd = design(d,'equiripple');
%     data=filter(Hd,data);


    %% or Bandpass filtering 3:
    [B_low,A_low] = butter(6,2*[fcl fch]./fs,'bandpass');
    data = filtfilt(B_low,A_low,data);
    %% Feature extracting
    [output_spectral_features,pxx,f,foct,spower,I,S] = spectral_features_3p(data,fs);
    OutputSheet(i+1,1)=cellstr(files(i).name);
    OutputSheet(i+1,2:14)=num2cell(output_spectral_features');
    OutputSheet(i+1,15:20)=num2cell(mfcc_mean);
end

writecell(OutputSheet,'OutputSheet_FaeBandpass_100_1000.xls');








fcl=125;  %% Low cutoff frequency
fch=625;    %% High cutoff frequency

OutputSheet=cell(size(files,1)+1,21);
OutputSheet(1,1:21)={'name','meanPSD','stdPSD','medPSD','bw','p25','p75','IQR','TP','p125_250','p250_375','p375_500','p500_625','spectrum_slope2','r_square2', 'mfcc1_mean', 'mfcc2_mean', 'mfcc3_mean', 'mfcc4_mean', 'mfcc5_mean', 'mfcc6_mean'};

%% Doing for each of the records ... 
for i=1:size(files,1)
    [data,fs]=audioread(['Samples\' files(i).name]);    %% read the record
    
    %% MFCC 
    mfcc_mean = mean(mfcc(data,fs, 'NumCoeffs', 6,'LogEnergy','Ignore'));

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
    %% Feature extracting
    [output_spectral_features,pxx,f,foct,spower,I,S] = spectral_features_4p(data,fs);
    OutputSheet(i+1,1)=cellstr(files(i).name);
    OutputSheet(i+1,2:15)=num2cell(output_spectral_features');
    OutputSheet(i+1,16:21)=num2cell(mfcc_mean);
end

writecell(OutputSheet,'OutputSheet_DavBandpass_125_625.xls');




OutputSheet=cell(size(files,1)+1,21);
OutputSheet(1,1:21)={'name','meanPSD','stdPSD','medPSD','bw','p25','p75','IQR','TP','p125_250','p250_375','p375_500','p500_625','spectrum_slope2','r_square2', 'mfcc1_mean', 'mfcc2_mean', 'mfcc3_mean', 'mfcc4_mean', 'mfcc5_mean', 'mfcc6_mean'};
%% Doing for each of the records ... 
for i=1:size(files,1)
    [data,fs]=audioread(['Samples\' files(i).name]);    %% read the record
    
    %% MFCC 
    mfcc_mean = mean(mfcc(data,fs, 'NumCoeffs', 6,'LogEnergy','Ignore'));

%     % Bandpass filtering 1 :     
%     dlp=fdesign.lowpass('Fp,Fst,Ap,Ast', .9*fch, 1.1*fch, 5, 60,fs);
%     Hdlp = design(dlp,'equiripple');
%     % fvtool(Hdlp)
%     % measure(Hd)
% 
%     dhp = fdesign.highpass('Fst,Fp,Ast,Ap', .9*fcl, 1.1*fcl, 60, 5, fs);
%     Hdhp = design(dhp,'equiripple');
%     % fvtool(Hdhp)
%     % measure(Hdhp)
% 
%     data=filter(Hdlp,(filter(Hdhp,data)));

    %% or Bandpass filtering 2 :
%     d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',.9*fcl,1.1*fcl,0.9*fch,1.1*fch,60,5,60,fs);
%     Hd = design(d,'equiripple');
%     data=filter(Hd,data);


    %% or Bandpass filtering 3:
    [B_low,A_low] = butter(6,2*[fcl fch]./fs,'bandpass');
    data = filtfilt(B_low,A_low,data);
    %% Feature extracting
    [output_spectral_features,pxx,f,foct,spower,I,S] = spectral_features_4p(data,fs);
    OutputSheet(i+1,1)=cellstr(files(i).name);
    OutputSheet(i+1,2:15)=num2cell(output_spectral_features');
    OutputSheet(i+1,16:21)=num2cell(mfcc_mean);
end

writecell(OutputSheet,'OutputSheet_FaeBandpass_125_625.xls');









