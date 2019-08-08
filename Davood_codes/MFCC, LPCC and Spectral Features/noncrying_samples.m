function [sample_num_minimum_noncrying_longest, sample_num_minimum_noncrying_longest_9sec, flag, y, fn, segmentTimeNonCry] = noncrying_samples(x,Fs)

%% INPUT and OUTPUT
% Gives the sample numbers of the minimum longest noncrying segments
% x = raw audio signal read as x
% Fs = sampling frequency
% sample_num_minimum_noncrying_longest = sample numbers of the minimum longest noncrying segments
% flag = longest noncrying segment > 4.5s => flag > 0 ok then flag = 0 no
% y = raw signal resampled
% fn = sampling frequency of y
% segmentTimeNonCry = time in second of non crying segment represented in tab


%% INITIALISATION

if (length(x)<60*Fs)
    tempLen = 60*Fs - length(x) + 1;
    x = [x;zeros(tempLen,1)]; % 0-padding to have 1 min minimum
end

x = x(1*Fs:60*Fs); % We take only the first minute of the signal (default choice)
y = resample(x,1,10); % Resampling of the signal by 10 (antialiasing)
fn = Fs/10; % New sampling frequency
dur = length(y)/fn;  % duration of the signal resampled


%% FREQUENCY ANALYSIS

% Producing the spectrogram and obtaining the frequencies, power and times
% of the audio sample

nfft = 1000; % fft points
window = ones(1,3*fn); % window = 3s
overlap = 0.5 * length(window); % Overlap = 1.5s

[s,f,t] = spectrogram(y, window, overlap, nfft, 'power', fn);
s_magnitude = abs(s); % magnitude of the power


%% POWER ANALYSIS

% Obtaining the power ratio for >350hz
f_greater_than_350hz = find(f>350); % to obtain only frequency > 350 Hz
sum_whole_column = sum(s_magnitude);

if (501<f_greater_than_350hz(end)) % 501 value obtained by Fatema's training
    sum_greater_than_350hz = sum(s_magnitude( f_greater_than_350hz(1)-1:501, 1:end ) );
else
    sum_greater_than_350hz = sum(s_magnitude(f_greater_than_350hz(1)-1:f_greater_than_350hz(end),1:end));
end

power_ratio = sum_greater_than_350hz./sum_whole_column;


%% DETERMINATION OF CRYING AND NON CRYING SEGMENT

% Obtaining the crying and non crying segments of the raw signal based on the power threshold

% Choice of the threshold determined with Fatema's training
index_crying = find(power_ratio(1,:)>0.5235);%power ratio when>0.5235 is crying, applying on raw signal
index_noncrying = find(power_ratio(1,:)<0.5235);%power ratio when<0.5235 is non-crying, applying on raw signal


% Creating a time vector (default choice : from 0s to 58.5s with overlap of
% 1.5s and window of 3s)
timeNumb = 0:1.5:floor(dur);
tempStr = strcat(num2str(timeNumb(1)),'-',num2str(timeNumb(1+2)));
timeString = {tempStr};

for ii=2:length(timeNumb)-2
    temp = strcat(num2str(timeNumb(ii)),'-',num2str(timeNumb(ii+2)));
    timeString = [timeString, {temp}];
end


% Making sure the length of the variables are even
if mod(length(power_ratio), 2) == 1 % x is odd
    power_ratio = power_ratio(:,1:(end-1)); % x is even
end
if mod(length(t), 2) == 1 % x is odd
    timeString = timeString(:,1:(end-1)); % x is even
end


%% DISPLAY OF CRYING AND NON CRYING SEGMENT

%producing merged cells containg the power ratios and corresponding time
%intervals
pr = power_ratio(1,:);
pr = num2cell(pr'); % converts numeric array to cell
C = cat(1,pr(:),timeString(:)); % power ratio then time
B = reshape(C,[],2); % power ratio and corresponding time intervals


%% CRYING

indexcell_crying = num2cell(index_crying); % index of crying segment eg : 1 2 4 5 6
filledCells = ~cellfun(@isempty,indexcell_crying); % calculating number of columns in the indexcell
columns_crying = sum(filledCells,2); % number of crying segment

% Produce the crying times
cry_times = num2cell(zeros(1,columns_crying)); % initialisation
for k = 1:columns_crying
    cry_times(k) = B(indexcell_crying{1,k},2);
end


%% NON CRYING

indexcell_noncrying = num2cell(index_noncrying); % index of non crying segment eg : 1 2 4 5 6

flag = 0; % to detect if a segment > 4.5 sec (by default no)
sample_num_minimum_noncrying_longest = 0; % initialization
sample_num_minimum_noncrying_longest_9sec = 0; % initialization
segmentTimeNonCry = [];

if (length(indexcell_noncrying)>1) % if there are non crying segment
    
    % Find consecutive segments (Mathieu)
    
    diffVector = [index_noncrying(1) diff(index_noncrying)]; % A : difference between cell n and n-1 of non crying segment
    mask = diffVector; % goal of mask : to detect consecutive segment (overlap) via 1
    % eg mask : -1 1 1 1 1 -1 5 6 -1 1 1 1 -1 => two consecutive segments
    % filled of 1 (-1 represents borders of segments)
    
    % detection of consecutive segments with mask
    if (diffVector(end) == 1)
        mask(end) = -1; % border consecutive segment : -1
    end
    if (diffVector(2) == 1)
        mask(1) = -1; % border consecutive segment : -1
    end
    for ii=2:length(diffVector)-1
        if (diffVector(ii) == 1 && diffVector(ii+1) == 1)
            mask(ii) = 0; % inside consecutive segment : 0
        else
            if (diffVector(ii) == 1 && (diffVector(ii+1) ~= 1)) || (diffVector(ii) ~= 1 && (diffVector(ii+1) == 1))
                mask(ii) = -1; % border consecutive segment : -1
            end
        end
    end
    
    indConsecutive = index_noncrying(find(mask==-1)); % find borders
    indOther = index_noncrying(find(mask>0)); % find non consecutive segment
    indConsecutiveTime = reshape((indConsecutive), 2, []); % reshaped index
    indConsecutiveTime(1,:) = (indConsecutiveTime(1,:) - 1).* 1.5; % start index in second (x1.5) (consecutive segment)
    indConsecutiveTime(2,:) = (indConsecutiveTime(2,:) - 1).* 1.5 +3; % end index in second (x1.5) (consecutive segment)
    indOtherTime = [(indOther-1).*1.5 ; (indOther-1).*1.5 + 3]; % start and end index in second (x1.5) (non consecutive segment)
    segmentTimeNonCry = [indConsecutiveTime indOtherTime]; % merge index
    segmentTimeNonCry = (sortrows(segmentTimeNonCry', 1,'ascend')); % list of whole time segment sorted
    % disp(segmentTimeNonCry);
    
    filledCells = ~cellfun(@isempty,indexcell_noncrying); % calculating number of columns in the indexcell
    columns_noncrying = sum(filledCells,2); % number of non crying segment
    
    
    % Produce the non crying times
    
    noncry_times = num2cell(zeros(1,columns_noncrying)); % initialisation
    for k = 1:columns_noncrying
        noncry_times(k) = B(indexcell_noncrying{1,k},2);
    end
%     disp('Non cry time');
%     disp(noncry_times);
    
    
    %%
    % Extracting the individual times
    out=cell2mat(cellfun(@str2num,strrep(noncry_times,'-',' '),'un',0)); % index decompose eg 1.5-3 3-4.5 => 1.5 3 3 4.5
    
    
    % Find consecutive segments (Mathieu)
    p2=find(diff(out)==-1.5);
    q2=[p2;p2+1];
    diff_2=out(q2);%-1.5s
    
    %     disp('diff2');
    %     disp(diff_2);
    %     disp('Segment');
    %     disp(segmentTimeNonCry');
    %     disp('Diff Segment');
    %     disp(diff(segmentTimeNonCry'));
    
    % Detect the longest segment
    [tempMax,ind] = max(diff(segmentTimeNonCry')); % difference between index in second for the start and the end of ech segment
    tempTime = segmentTimeNonCry(ind,1):1.5:segmentTimeNonCry(ind,2); % new temporal axis adapted to the longest segment
    %     disp('ind');
    %     disp(ind);
    %     disp('1');
    %     disp(segmentTimeNonCry(ind,1));
    %     disp('2');
    %     disp(segmentTimeNonCry(ind,2));
    
    
    if (size(diff_2,1) == 2) % if there are consecutive segments
        flag = 1; % flag to know there are a non crying segment
        
        tempIndexStart = find(diff_2(2,:) == segmentTimeNonCry(ind,1)+1.5); % find the beginning of the longest sample
        tempIndexEnd = find(diff_2(1,:) == segmentTimeNonCry(ind,2)-1.5); % find the end of the longest sample
        
        
        sample_num_minimum_noncrying_longest = diff_2(2*tempIndexStart:2*tempIndexEnd-1)*fn; % row (thus x2 for the index) coefficients (in sample)
        if isempty(sample_num_minimum_noncrying_longest)
            flag = 0;
        end
        
        if (tempMax-3 >=4.5) % if longest segment > 4.5 sec
            
            flag = 2; % flag to know there are a segment > 4.5 sec
            
            middleTime = length(tempTime)/2;
            startTime = max([0, ceil(middleTime)-1]);
            disp(startTime);
            finishTime = min([length(tempTime), startTime + 3]);
            tempIndexStart9sec = find(diff_2(2,:) == tempTime(startTime));
            tempIndexEnd9sec = find(diff_2(1,:) == tempTime(finishTime));
            
%             disp('diff');
%             disp(diff_2);
%             disp('tempSatrt');
%             disp(tempTime);
%             disp(tempTime(startTime));
%             disp(tempTime(finishTime));
%             disp('index');
%             disp(tempIndexStart9sec);
%             disp(tempIndexEnd9sec);
            
            
            minimum_seg_longest_noncrying_longest_9sec = y(tempTime(startTime)*fn+1:1:tempTime(finishTime)*fn); % 4.5-sec segment takes in the middle of longest segment
            audiowrite('minimum_seg_longest_noncrying_9sec.wav',minimum_seg_longest_noncrying_longest_9sec,fn); % writting of this segment
            sample_num_minimum_noncrying_longest_9sec = diff_2(2*tempIndexStart9sec:2*tempIndexEnd9sec-1); % row (thus x2 for the index) coefficients (in sample)
%                 disp('EXAMPLE');
%                 disp('MIDDLE OF A SEGMENT:')
%                 disp([tempTime(startTime) tempTime(finishTime)]);
%                 disp(tempTime(finishTime)-tempTime(startTime));
        end
        
        %% if one wishes to find out the clipped segment
        % Otherwise comment out this section
        % Minimum longest noncrying continuous
        % Unfiltered
        
        % disp(length(tempTime));
        minimum_seg_longest_noncrying = y(tempTime(1)*fn+1:1:tempTime(end)*fn); % take the whole longest segment
        audiowrite('minimum_seg_longest_noncrying.wav',minimum_seg_longest_noncrying,fn); % writting of this segment
        % disp('WHOLE SEGMENT:')
        % disp([tempTime(1) tempTime(end)]);
        % disp(tempTime(end)-tempTime(1));
        % pause(2);
    end
end


end

