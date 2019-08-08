function [sample_num_minimum_noncrying_longest, y, fn] = noncrying_samples2(x,Fs)
%%
% Gives the sample numbers of the minimum longest noncrying segments
% input = raw audio signal read as x
% Fs = sampling frequency
% output = sample numbers of the minimum longest noncrying segments

%%
close all

x = x(Fs*1:Fs*60);
y = resample(x,1,2); fn = Fs/2; t = (0:1/fn)';

  dur = length(y)/fn;%duration of the signal

  %%
%producing the spectrogram and obtaining the crequencies, power and times
%of the audio sample
% [s,f,t] = spectrogram(x,window,noverlap,f,fs) 
% returns the spectrogram at the cyclical frequencies specified in f.
%window = 3s, Overlap = 1.5s, fft points = 1000
[s,f,t] = spectrogram(y,66150,33075,1000,'power', fn);

 s_magnitude = abs(s);%magnitude of the power

 %%
 %obtaining the power ratio for 350hz
 sum_whole_column = sum(s_magnitude);
 sum_greater_than_350hz = sum(s_magnitude( 17:501, 1:end ) );
power_ratio = sum_greater_than_350hz./sum_whole_column;

%%
 %obtaining the crying and non crying segments of the raw signal based on the power threshold
 index_crying = find(power_ratio(1,:)>0.5235);%power ratio when>0.5235 is crying, applying on raw signal

 index_noncrying = find(power_ratio(1,:)<0.5235);%power ratio when<0.5235 is non-crying, applying on raw signal
 
 % total number of samples in the raw signal
total_samples = [0,dur*fn];

%creating a time vector
    time =  {'0-3',
        '1.5-4.5',
        '3-6',
        '4.5-7.5',
        '6-9',
        '7.5-10.5',
        '9-12',
        '10.5-13.5',
        '12-15',	
        '13.5-16.5',
        '15-18',	
        '16.5-19.5',
        '18-21',	
        '19.5-22.5',
        '21-24',	
        '22.5-25.5',
        '24-27',	
        '25.5-28.5',	
        '27-30',	
        '28.5-31.5',
        '30-33',	
        '31.5-34.5',
        '33-36',	
        '34.5-37.5',
        '36-39',	
        '37.5-40.5',	
        '39-42',	
        '40.5-43.5',	
        '42-45',	
        '43.5-46.5',
        '45-48',	
        '46.5-49.5',
        '48-51', 
        '49.5-52.5',
        '51-54',	
        '52.5-55.5',	
        '54-57',
	'55.5-58.5'};

%making sure the length of the variables are even
if mod(length(power_ratio), 2) == 1 % x is odd
      power_ratio = power_ratio(:,1:(end-1));
else
 power_ratio = power_ratio; % x is even
end
if mod(length(t), 2) == 1 % x is odd
      t = t(:,1:(end-1));
else
 t = t; % x is even
end

%%
%producing merged cells containg the power ratios and corresponding time
%intervals
    pr = power_ratio(1,:);
    pr = num2cell(pr');%converts numeric array to cell
    C = cat(1,pr(:),time(:));
    B = reshape(C,[],2); %power ratio and corresponding time intervals

    %%
% %for crying

indexcell_crying = num2cell(index_crying);
filledCells = ~cellfun(@isempty,indexcell_crying);%calculating number of columns in the indexcell
columns_crying = sum(filledCells,2);

% produce the crying times

for k = 1:(columns_crying)  
cry_times(k) = B(indexcell_crying{1,k},2);
% disp(cry_times(k));
end
%%
% %for non-crying
indexcell_noncrying = num2cell(index_noncrying);
filledCells = ~cellfun(@isempty,indexcell_noncrying);%calculating number of columns in the indexcell
columns_noncrying = sum(filledCells,2);
% produce the non crying times
for k = 1:columns_noncrying    
noncry_times(k) = B(indexcell_noncrying{1,k},2);
disp('non cry');
disp(noncry_times(k));
end
%%
%extracting the individual times
out=cell2mat(cellfun(@str2num,strrep(noncry_times,'-',' '),'un',0));
columns_noncry_times = length(out);
%%
%sample numbers of the non crying segments
for i = 1:columns_noncry_times
    for j = 1:2:columns_noncry_times
sample_begin_noncrying(i,j) = out(1,j)*fn;
sample_end_noncrying(i,j) = out(1,j+1)*fn;
    end
    sample_num_noncrying = [reshape(sample_begin_noncrying,[],1) reshape(sample_end_noncrying,[],1)];
    
end

%%
%sample numbers for continuous 9s
%taking segments with overlapping
difference_times_noncry1 =diff(out);

p2=find(diff(out)==-1.5);
q2=[p2;p2+1];
diff_2=out(q2);%-1.5s

%%
%final output
%read as column2-column1; column4-column3 and so on
disp(length(diff_2));
disp(diff_2);
sample_num_minimum_noncrying_longest = diff_2(3:14)*fn;%taking the middle segment
%% if one wishes to find out the clipped segment
% otherwise comment out this section
 %minimum longest noncrying continuous
 %unfiltered
 disp(sample_num_minimum_noncrying_longest(1,2)/fn);
 disp(sample_num_minimum_noncrying_longest(1,end-1)/fn);
 minimum_seg_longest_noncrying = x(sample_num_minimum_noncrying_longest(1,2) : sample_num_minimum_noncrying_longest(1,end-1));
%  plot(minimum_seg_longest_noncrying);
 audiowrite('minimum_seg_longest_noncrying.wav',minimum_seg_longest_noncrying,fn);
 
 sample_num_minimum_noncrying_longest_9sec = diff_2(9:14)*fn;%taking the middle segment
 minimum_seg_longest_noncrying_9sec = x(sample_num_minimum_noncrying_longest_9sec(1,2) : sample_num_minimum_noncrying_longest_9sec(1,end-1));
%  plot(minimum_seg_longest_noncrying);
 audiowrite('minimum_seg_longest_noncrying_9sec.wav',minimum_seg_longest_noncrying_9sec,fn);
end

% %%
% %final output
% %read as column2-column1; column4-column3 and so on
% sample_num_minimum_noncrying_longest = diff_2*fn;
% if length(diff_2)>=9
%     middle = length(diff_2)/2;
%     start = ceil(middle-4);
%     finish = start + 8;
%     sample_num_minimum_noncrying_longest_9sec = diff_2(start:finish)*fn;%taking the middle segment
%     minimum_seg_longest_noncrying_longest_9sec = x(sample_num_minimum_noncrying_longest_9sec(1,2) : sample_num_minimum_noncrying_longest_9sec(1,end-1));
%     %  plot(minimum_seg_longest_noncrying);
%     audiowrite('minimum_seg_longest_noncrying_9sec.wav',minimum_seg_longest_noncrying_longest_9sec,fn);
% end
% %% if one wishes to find out the clipped segment
% % otherwise comment out this section
%  %minimum longest noncrying continuous
%  %unfiltered
%  minimum_seg_longest_noncrying = x(sample_num_minimum_noncrying_longest(1,2) : sample_num_minimum_noncrying_longest(1,end-1));
% %  plot(minimum_seg_longest_noncrying);
%  audiowrite('minimum_seg_longest_noncrying.wav',minimum_seg_longest_noncrying,fn);
% end
