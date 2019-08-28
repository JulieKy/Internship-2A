clear all;

truc=9;
x=linspace(1, truc, truc);
N=length(x);

window=3;
overlap=1;

%% SECTIONS
% Initialisation
pos=1;
section_set=[];

% All sections except the last one 
while pos<=(N-window+1)
    section=x(pos:pos+window-1)';
    pos=pos+window-overlap;
    section_set=[section_set, section];
end

% Last section
last_sec=x(pos:end); % Can have a different size

%% POWERBAND
label_section=mean(section_set)>3; % For all sections except the last one
label_section_last=mean(last_sec)>3; % Last section


%% REMETTRE p AVEC MEME NOMBRE QUE X AVEC OVERLAP

% -- First section
first_section=x(1:window-overlap);
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
