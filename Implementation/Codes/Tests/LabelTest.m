% Labels the recordings and give a matrix with the 3 observators (dire mean
% et tout ca): output=labels & input=window et overlap samples, observators

%% INITIALISATION

% -- Data
observators=2;
samples=37;

% -- Parameters
window=3;
overlap=25/100;   % A METTRE DANS APPEL DE FONCTION

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
