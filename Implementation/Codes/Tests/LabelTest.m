% Labels the recordings and give a matrix with the 3 observators (dire mean
% et tout ca): output=labels & input=window et overlap samples, observators
clear all, close all, clc, dbstop if error;
addpath(genpath('..\..\')); % to have access to sample folder
path = pwd; % current path
data_dir=[path,'\..\..\Data\Samples_Belle\'];
dinfo = dir(data_dir);
names_cell1 = {dinfo.name};
% choose a valid file name
n_section=0;
for i=1:size(names_cell1,2)
    if length(names_cell1{i})>2
        n_section=n_section+1;
        names_cell{n_section}=names_cell1{i};
    end
end

close all; % to close previous figures opened in computation of spectral features
tempName=names_cell{1};
% to follow which file is red
disp('READ');
disp(tempName);

[x,Fs]= audioread([path,'\..\..\Data\Samples_Belle\',tempName]); % read current file

%% resampling to 4000 Hz
xs=resample(x,4000,Fs);
fn=4000;

%% ------------------- Shorten the signals to 60s A METTRE DANS LE CODE ---------------------
time_sample=60;
xss=xs(1:time_sample*fn,1);

%% -------------------------------------

%% INITIALISATION

% -- Data
observators=2;
samples=37;

% -- Parameters
window=1;
overlap=0;   % A METTRE DANS APPEL DE FONCTION

% -- Variables
end_sample=60; % End of the signal (hypotesis: length of the signal=60s)


%% LABELLING THE BANK OF SAMPLES
for samp=1:samples % Number of samples
    for obs=1:observators % Number of observators
        
        
        start=0; % Departure of the section
        end_w=start+window; % End of the window section
        label_signal=[]; % Labels of one signal
        
        
        %% READING A TEXT FILE
        % The files are named as follows: a_b.txt with a=observator and b=sample
        
        % -- Openning the file
        file_name=sprintf('Labels\\%d_%d.txt', obs, samp);
        fid=fopen(file_name);
        
        % -- Extract the data
        fmt=['%n', '%n', '%s'];
        file=textscan(fid,fmt);
        
        % -- Close the file
        fclose(fid);
        
        % -- Reading the data
        file_time_s=cell2mat(file(1,1)); % The start time of a section CS or NCS (colon 1 in the text file)
        file_time_e=cell2mat(file(1,2)); % The ending time of a section CS or NCS (colon 2 in the text file)
        file_label=file{3}; % The label of a section CS or NCS (colon 3 in the text file)
        
        
        %% LABELLING
        l=1; % Line of the text file
        
        % -- For an entire signal
        while start<end_sample-window+window*overlap
            
            flag=0; % To know if there are different sections (CS&NCS) in one window
            label_window=[]; % The label in a window, with a precision of 10^-2
            
            time_start=file_time_s(l); % The start time of a section
            time_end=file_time_e(l); % The ending time of a section
            
            % -- For one window
            while end_w>time_end
                flag=1; % There are different sections (CS&NCS) in one window
                label_section=file_label{l}; % The label during a section (colon 3 in the text file)
                
                % -- Filing label_window with a number of 1(CS) or 0(NCS) corresponding to the section duration
                delta_t=round((time_end-start)*100); % The duration of a section type in a window rounded to 10^-2, *100
                if strcmp(label_section, 'CS')==1 % Crying Section: 1
                    label_window=[label_window ones(1, delta_t)];
                elseif strcmp(label_section, 'NCS')==1 % Non Crying Section: 0
                    label_window=[label_window zeros(1, delta_t)];
                else
                    error('The label must be ''CS'' or ''NCS''');
                end
                
                % -- Next section
                l=l+1;
                time_end=file_time_e(l);
            end
            
            % -- Ending the window (when end_w<time_end but end_w>start)
            label_section=file_label{l};
            
            % Different (flag=1) or no different (flag=0) section types in one window
            if flag==1
                delta_t=round((end_w-file_time_e(l-1))*100);
            else
                delta_t=round((end_w-start)*100);
            end
            
            % Filing label_window with a number of 1(CS) or 0(NCS) corresponding to the section duration
            if strcmp(label_section, 'CS')==1
                label_window=[label_window ones(1, delta_t)];
            elseif strcmp(label_section, 'NCS')==1
                label_window=[label_window zeros(1, delta_t)];
            else
                error('The label must be ''CS'' or ''NCS''');
            end
            
            % -- Filling the vector of signal labels by establishing the window section as 'CS' or 'NCS'
            label_signal=[label_signal mean(label_window)>0.5]; % CS=1 and NCS=0
            
            % -- New window section
            start=start+window-window*overlap;
            end_w=start+window;
        end
        
        % -- Resulting labels of the bank of signals
        labels(obs, :, samp)=label_signal;
    end
end

%% KAPPA CALCULATION
%  -- Initialisation
coef_KAPPA=zeros(1, samples);

label_final=[]; % A METTRE EN HAUT!!

for samp=1:samples % Number of samples
    
    %  -- Counting the number of raters who agreed on a category (CS and NCS)
    x=labels(:, :, samp); % Matrix shape: (labels of one sample, observators)
    nb_CS=sum(x); % Number of raters who agreed on CS (vector)
    nb_NCS=observators-nb_CS; % Number of raters who agreed on NCS (vector)
    y=[nb_CS', nb_NCS'];
    
    % -- Calculation for each sample
    coef_KAPPA(1,samp)=fleiss(y);
    
    %% LABELLING BY FINDING A CONCENSUS WHEN RATERS DON'T AGREE
    % -- Determining one label per window, putting the 'CS' label when there are at least 2/3 of agreement
    threshold_agreement=2/3;
    label_final=[label_final; (nb_CS/observators)>threshold_agreement];
end

%% POWER RATIO

i=1; % NUMERO DE SAMPLE A DONNER!!

N = length(xss);
time_axis = (1:N)/fn;


% PEUT ETRE A METTRE DANS FONCTION
% -- Finding the location of 'CS' and 'NCS'
CS_locs=find(label_final(i,:)==1);
NCS_locs=find(label_final(i,:)==0);

% -- Start time of the labels (for each window)
time_CS_start=CS_locs*(window-overlap);
time_NCS_start=NCS_locs*(window-overlap);

% -- Start sample of the labels (for each window)
sample_CS_start=time_CS_start*fn;
sample_NCS_start=time_NCS_start*fn;
label_duration=window*fn; % Number of samples in a window

pass_band=[0:1000]; % PAS ICI
f_band_width=100;  % Width of frequency bands
% -- Periodogram for each CS sections
pr_CS=zeros(length(CS_locs),floor(pass_band(end)/f_band_width)); % Power ratio
for n_section=1:length(CS_locs)
    CS_section=xss(sample_CS_start(n_section):sample_CS_start(n_section)+label_duration);
    % FAIRE UN PLOT
    [pxx,f] = Welch_periodogram(CS_section, fn, pass_band);  % A MIEUX FAIRE
    pxx_sum=sum(pxx);
    for h = pass_band(1)+1 : floor(pass_band(end)/f_band_width) % A MIEUX FAIRE
        a=1; i=1;
        while f(i)<100*h
            i=i+1;
        end
        pr_CS(n_section,h)=mean(pxx(a:i-1))/pxx_sum;
        a=i;
    end
end

pr_CS_mean=mean(pr_CS);

% -- Periodogram for each NCS sections
pr_NCS=zeros(length(NCS_locs),floor(pass_band(end)/f_band_width)); % Power ratio

% For each NCS section
for n_section=1:length(NCS_locs)
    NCS_section=xss(sample_NCS_start:sample_NCS_start+label_duration); % Section NCS
    [pxx,f] = Welch_periodogram(NCS_section, fn, pass_band);  % Welch Periodogram of the section
    f_interval=pass_band(end)/length(f); % Frequency interval in the f vector
    band_interval=round(f_band_width/f_interval)+1; % Index interval for a frequency band in the f vector
    pxx_sum=sum(pxx); % Sum of the pxx
    
    % For each frequency band
    for n_band = pass_band(1) : floor(pass_band(end)/f_band_width)-1
        % pr_NCS(n_section,n_band+1)=mean(pxx(band_interval*n_band+1:band_interval*(n_band+1)))/pxx_sum; % Power ratio
        % A VOIR !!
    end  
end

pr_NCS_mean=mean(pr_NCS); % Mean of the power ratio for all NCS section of a sample