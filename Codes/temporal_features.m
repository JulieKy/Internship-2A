% This function computes temporal features 

function [output_temporal_features] = temporal_features(x,fn)
%% INPUT AND OUTPUT

% -- INPUTS
% x 'audio signal'
% fn 'sampling frequency

% -- OUTPUTS =  temporal features including

%% VARIABLES
N = length(x);
time_axis = (1:N)/fn;

%% ZERO CROSSING RATE
ZCR=sum(abs(diff(sign(x))/2))/length(x);

%% DISPLAY
 figure,
 plot(time_axis, x);
 title('Temporal representation');
 
%% OUTPUTS   
  output_temporal_features=ZCR;

end

