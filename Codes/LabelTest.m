% Labels the recordings and give a matrix with the 3 observators (dire mean
% et tout ca)

%% INITIALISATION

% -- Data
observators=1; % a mettre dans l'appel de focntion
samples=1; %idem

% -- Parameters
window=3;
overlap=25/100;
Label_w=[];

% -- Sections
d=0;            % Departure of the section
e=d+window;     % End of the section
E=60;           % End of the signal (hypotesis: length of the signal=60s)


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
        while d<E-window+window*overlap % For an entire signal
            time_end=section_time(l); % The ending time of a section CS or NCS (colon 2 in the text file)
            while e>time_end % One window section
                label_w=section_label{l};
                if label_w=='CS' % Crying Section
                    Label_w=[Label_w ones(20,1)];
                elseif label_w=='NCS' % Non Crying Section
                    Label_w=[Label_w zeros(20,1)];
                else 
                    error('The label must be ''CS'' or ''NCS''');
                end
                    
            end
            d=d+window-window*overlap;
        end
        
    end
end

