function [ X ] = read_samples( folder_path, time_sample, fn )
%read_samples: Read the samples of a folder and return all the folder in a matrix. Each line number corresponds to the sample number.

%% INPUTS AND OUTPUTS
% --- Inputs ---
% folder_path: Absolute path
% --- Outputs ---
% X: Matrix with all the samples. Each line number corresponds to the sample number.

%% -- Samples' names initialisation
dinfo = dir(folder_path);
names_cell1 = {dinfo.name};

% Choose a valid file name
j=0;
for i=1:size(names_cell1,2)
    if length(names_cell1{i})>2
        j=j+1;
        names_cell{j}=names_cell1{i};
    end
end
lengthTot=j;

%% -- Reading 
X=zeros(lengthTot, time_sample*fn); % All the signals stored in this matrix

for i = 1:lengthTot % loop to have all recording
    
    close all; % Close previous figures opened in computation
    
    % Name of the sample
    tempName=names_cell{i};
    disp('READ');
    disp(tempName);
    
    % Get the number of the recording by removing the '.mp3'
    strMP3 = sprintf('%s',tempName);
    ind=strfind(strMP3,'.');
    signal_n = str2num(strMP3(1:ind-1));
    if (signal_n~=34 && signal_n~=48 && signal_n~=49)
    
    % Reading the sample
    sample_path=sprintf('%s\\%s', folder_path, tempName);
    [x,Fs]= audioread(sample_path); % read current file
    
    % Resampling to 4000 Hz
    xs=resample(x,fn,Fs);
    
    % Shorten the signals to 60s
    xss=xs(1:time_sample*fn,1);
    
    % Fill X with all samples
    X(signal_n, :)=xss;
    end
end

end
