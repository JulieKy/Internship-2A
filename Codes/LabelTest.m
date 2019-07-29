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
Label_w=[];
Label_w_mean=[];

% -- Sections
d=0; % Departure of the section
end_w=d+window; % End of the window section
end_sample=60; % End of the signal (hypotesis: length of the signal=60s)


%% READING A TEXT FILE
% The files are named as follows: a_b.txt with a=observator and b=sample

% -- Openning the file
for o=1:observators % Number of observators
    for s=1:samples % Number of samples
        file_name=sprintf('Labels\\%d_%d.txt', o, s);
        fid=fopen(file_name);
        
        % -- Extract the data
        fmt=['%n', '%n', '%s'];
        file=textscan(fid,fmt);
        
        % -- Close the file
        fclose(fid);
        
        % -- Reading the data
        section_time=cell2mat(file(1,2));
        section_label=file{3};
        
        
        %% LABELLING
        l=1; % Line of the text file
        
        % -- For an entire signal
        while d<end_sample-window+window*overlap
            time_end=section_time(l); % The ending time of a section CS or NCS (colon 2 in the text file)
                        
            % -- For one window
            while end_w>time_end
                label_w=section_label{l}; % The label during a section CS or NCS (colon 3 in the text file)
                
                % -- Filing Label_w with a number of 1(CS) or 0(NCS) corresponding to the section duration
                delta_t=round((section_time(l+1)-section_time(l))*100); % The duration of a section rounded to 10^-2, *100
                if strcmp(label_w, 'CS')==1 % Crying Section: 1
                    Label_w=[Label_w ones(1, delta_t)];
                elseif strcmp(label_w, 'NCS')==1 % Non Crying Section: 0
                    Label_w=[Label_w zeros(1, delta_t)];
                else
                    error('The label must be ''CS'' or ''NCS''');
                end
                
                l=l+1;
                time_end=section_time(l);
            end
            
            % -- Ending the window (when end_w<time_end but end_w>d)
            label_w=section_label{l};
            delta_t=round((section_time(l)-end_w)*100); % The duration of a section rounded to 10^-2, *100
            if strcmp(label_w, 'CS')==1
                Label_w=[Label_w ones(1, delta_t)];
            elseif strcmp(label_w, 'NCS')==1
                Label_w=[Label_w zeros(1, delta_t)];
            else
                error('The label must be ''CS'' or ''NCS''');
            end
            
            % -- Establishing the window section as 'CS' or 'NCS'
            Label_w_mean=[Label_w_mean mean(Label_w)>0.5]; % CS=1 and NCS=0
            
            % -- New window section
            d=d+window-window*overlap;
            end_w=d+window;
        end
    end
end

