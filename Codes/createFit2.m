function [fitresult, gof] = createFit2(f, power)
%CREATEFIT(F,POWER)
%  Create a fit.
%%producing linear regression line fit to logarithmic scale psd
%needs to be called by spectral_features function
%  Data for 'untitled fit 1' fit:
%      X Input : f
%      Y Output: power
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( f, power );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );



