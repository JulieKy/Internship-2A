% This function computes temporal features

function [output_temporal_features] = temporal_features(x,fn, tempName)
%% INPUT AND OUTPUT

% -- INPUTS
% x 'audio signal'
% fn 'sampling frequency

% -- OUTPUTS =  temporal features including
% ZRC: Zero Crossing Rate

%% VARIABLES
N = length(x);
time_axis = (1:N)/fn;

%% ZERO CROSSING RATE
ZCR=sum(abs(diff(sign(x))/2))/length(x);

%% DISPLAY
fig=figure;
plot(time_axis, x);

% Get the number of the recording by removing the '.mp3'
strMP3 = sprintf('%s',tempName);
ind=strfind(strMP3,'.');
strNum = strMP3(1:ind-1);

% Title and legend
strTitle=sprintf('Recording %s',strNum);
title(strTitle,'fontsize',14,'interpreter','latex');
xlabel('Time (s)'),ylabel('Amplitude');

%% SAVE 
% Create a folder 
path = pwd;
pathFigTemporal = strcat(path, '\Figures\Time_figures');
if ~(exist(pathFigTemporal)) % test to create excel file or no
    disp('Creation temporal representation folder')
    mkdir Time_figures
end

% Save the figure
saveas(fig, fullfile(pathFigTemporal, strNum), 'png');

%% OUTPUTS
output_temporal_features=ZCR;

end

