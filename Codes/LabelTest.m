% Labels the recordings and give a matrix with the 3 observators (dire mean
% et tout ca)

%% INITIALISATION

% -- Data
observators=1; % a mettre dans l'appel de focntion
samples=1; %idem

% -- Parameters
% window=3;
% overlap=25/100;
window=1;
overlap=0;

label_signal=[]; % A METTRE POUR CHAQUE SIGNQUX

% -- Sections
start=0; % Departure of the section
end_w=start+window; % End of the window section
end_sample=60; % End of the signal (hypotesis: length of the signal=60s)


%% READING A TEXT FILE
% The files are named as follows: a_b.txt with a=observator and b=sample

% -- Openning the file
% for obs=1:observators % Number of observators
%     for samp=1:samples % Number of samples

obs=1;
samp=1;

file_name=sprintf('Labels\\%d_%d.txt', obs, samp);
fid=fopen(file_name);
strtxt=sprintf('READING %d_%d.txt', obs, samp);
disp(strtxt)

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
i=0;

% -- For an entire signal
while start<end_sample-window+window*overlap
    
    flag=0;
    
    i=i+1;
    strtruc=sprintf('--- Boucle signal: %d', i);
    disp(strtruc)
    
    label_window=[];
    
    time_start=file_time_s(l); % The start time of a section
    time_end=file_time_e(l); % The ending time of a section
    
    strts=sprintf('time_start1: %d', time_start);
    disp(strts)
    strtsaa=sprintf('start: %d', start);
    disp(strtsaa)
    strte=sprintf('time_end1: %d', time_end);
    disp(strte)
    stre=sprintf('end_w1: %d', end_w);
    disp(stre)
    
    % -- For one window
    while end_w>time_end
        flag=1;
        label_section=file_label{l}; % The label during a section CS or NCS (colon 3 in the text file)
        
        % -- Filing Label_w with a number of 1(CS) or 0(NCS) corresponding to the section duration
        delta_t=round((time_end-start)*100); % The duration of a section rounded to 10^-2, *100
        
        
        strdt1=sprintf('delta t window: %d', delta_t);
        disp(strdt1)
        
        if strcmp(label_section, 'CS')==1 % Crying Section: 1
            label_window=[label_window ones(1, delta_t)];
        elseif strcmp(label_section, 'NCS')==1 % Non Crying Section: 0
            label_window=[label_window zeros(1, delta_t)];
        else
            error('The label must be ''CS'' or ''NCS''');
        end
        
        strw=sprintf(' --- Boucle window %d', l);
        disp(strw)
        
        l=l+1;
        time_end=file_time_e(l);
        
        strts=sprintf('time_start2: %d', time_start);
        disp(strts)
        strte=sprintf('time_end2: %d', time_end);
        disp(strte)
        stre=sprintf('end_w2: %d', end_w);
        disp(stre)
    end
    
    % -- Ending the window (when end_w<time_end but end_w>d)
    
    strfw=sprintf('--- Fin window %d', l-1);
    disp(strfw)
    
    label_section=file_label{l};
    
    if flag==1
        delta_t=round((end_w-file_time_e(l-1))*100); % The duration of a section rounded to 10^-2, *100
    else
        delta_t=round((end_w-start)*100);
    end
    
    strdt2=sprintf('delta t fin: %d', delta_t);
    disp(strdt2)
    
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
%     end
% end

