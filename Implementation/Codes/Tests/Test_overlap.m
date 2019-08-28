clear all;

folder_path='C:\Users\julie\OneDrive\Documents\T2\Stage 2A\Stage_2A_Matlab\Implementation\Codes\..\Data\Database\';
time_sample=60;
fn=4000;
threshold=0.0026;
band=[296, 407];
window_training=3;
overlap_training=1;

window=window_training*fn;
overlap=overlap_training*fn;

%% READING FILES IN THE DATABASE

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

%% -- Initialisation storage matrix
X_ncs=zeros(lengthTot, time_sample*fn);

%% -- Reading
for i=1:lengthTot
    
    % Name of the sample
    tempName=names_cell{i};
    disp('READ - Database');
    disp(tempName);
    
    % Get the number of the recording by removing the '.mp3'
    strMP3 = sprintf('%s',tempName);
    ind=strfind(strMP3,'.');
    signal_n = str2num(strMP3(1:ind-1));
    
    % Reading the sample
    sample_path=sprintf('%s\\%s', folder_path, tempName);
    [x,Fs]= audioread(sample_path); % read current file
    
    % Resampling to 4000 Hz
    xs=resample(x,fn,Fs);
    
    % Shorten the signals to 60s if longer
    if length(xs)>time_sample*fn
        xss=xs(1:time_sample*fn,1);
    end
    
    
    %% SECTIONS
    % Initialisation
    pos=1;
    section_set=[];
    N=length(xss); 
    
    % All sections except the last one
    while pos<=(N-window+1)
        section=xss(pos:pos+window-1);
        pos=pos+window-overlap;
        section_set=[section_set, section];
    end
    
    % Last section
    last_sec=xss(pos:end); % Can have a different size
    
    %% POWERBAND
    label_section=bandpower(section_set, fn, band); % For all sections except the last one
    label_section_last=bandpower(last_sec, fn, band); % Last section
    
    
    %% REMETTRE p AVEC MEME NOMBRE QUE X AVEC OVERLAP
    
    % -- First section
    first_section=xss(1:window-overlap);
    labels(1:length(first_section))=repelem(label_section(1),length(first_section));
    pos=length(first_section)+1;
    
    % -- Other windows
    for i=2:length(label_section)
        l1=label_section(i-1);
        l2=label_section(i);
        
        % Overlap with priority to CS
        if (l1==1 || l2==1)
            labels(pos: pos+overlap-1)=ones(1,overlap);
        else
            labels(pos: pos+overlap-1)=zeros(1,overlap);
        end
        pos=pos+overlap;
        
        % Window without overlap
        labels(pos:pos+window-2*overlap-1)=repelem(l2,window-2*overlap);
        pos=pos+window-2*overlap;
    end
    
    % -- Last window
    
    % Last overlap
    if (label_section(end)==1 || label_section_last==1)
        labels(pos: pos+overlap-1)=ones(1,overlap);
    else
        labels(pos: pos+overlap-1)=zeros(1,overlap);
    end
    pos=pos+overlap;
    
    % Last section without overlap
    if length(last_sec)~=1 % If length=1, it means that there is only overlap (already taken into account)
        labels(pos: pos+length(last_sec)-overlap-1) = repelem(label_section_last, length(last_sec)-overlap);
    end
    
    
    %% CS REMOVING
    NCS=1-labels;
    xsc=xss((xss.*NCS')~=0); % Signal without CS
    
    %% STORAGE
    X_ncs(signal_n, 1:length(xsc))=xsc; % If CS were removed, X_ncs contains 0.
    
    %% LENGTH
    length_xsc(signal_n)=length(xsc);
    
    %% LABEL TRAINING
    %label_training(signal_n, :)=repelem(CS, round(window_training/window_annotated)); % Useful in 'display_CS_NCS_final.m'
    
end


%% SHORTEN SAMPLES WITH MINIMUM LENGTH
% Find the minimum length
min_length=min(length_xsc);
part1=floor(min_length/2);
part2=min_length-part1;

length_time=length_xsc./fn;
min_ok=10;
figure,
hax=axes;
x_axe=get(hax,'XLim');
plot(1:length(length_time), length_time, '*'); hold on
plot(1:length(length_time),ones(1,length(length_time))*10 , '--', 'Color', 'r');
hold off
title('Duration of Samples after Cry Removal')
xlabel('Samples')
ylabel('Duration [s]')


% Statistical study of lengths
length_xsc_time=length_xsc*60/length(xss);
mean_length=mean(length_xsc_time);
median_length=median(length_xsc_time);
p25_length=prctile(length_xsc_time,25);
p75_length=prctile(length_xsc_time,75);


% Shorten for each signal
X_ncs_final=zeros(lengthTot, min_length);
for signal_n=1:size(X_ncs,1)
    xsc1=X_ncs(signal_n, 1:length_xsc(signal_n)); % Signal with only NCS
    if length(xsc1)>min_length
        mid=floor(length(xsc1)/2);
        xsc_shorten=xsc1(mid-part1:mid+part2-1);
    else
        xsc_shorten=xsc1;
    end
    X_ncs_final(signal_n, :)=xsc_shorten;
end