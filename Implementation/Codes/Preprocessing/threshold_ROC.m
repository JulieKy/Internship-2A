function [ fpr, tpr, final_threshold ] = threshold_ROC( nb_thresholds, powerband )
%threshold_ROC: Based on powerband. Compute the ROC.

%% INPUTS AND OUTPUTS
%  -- Inputs --
% thresholds_n: Number of thresholds for the ROC
% powerband: Powerband of all sections, with the tag_section (0=NCS, 1=CS). Shape of the matrix: [tag_section, powerband].
% -- Outputs --
% thershold_final: The better threshold to determine CS and NCS

%% INITIALISATION
thresholds_test=linspace(0, max(powerband(:,2)), nb_thresholds); % Thresholds for the ROC

tpr=[]; % True Positive Ratio
fpr=[];  % False Positive Ratio

%% ROC COMPUTATION
% For each threshold
for i=1:nb_thresholds
    threshold=thresholds_test(i);
    
    % -- Compare powerband of each section to the threshold
    label_obtained=(powerband(:,2)>threshold)*2; % Hypothesis: 2 for CS when powerband>threshold; 0 for NCS when powerband<threshold
    
   % -- Compare the labels obtained to the real labels
   res=powerband(:,1) + label_obtained; % With 'real'/'obtained', 3=CS/CS; 1=CS/NCS; 2=NCS/CS; 0=NCS/NCS; 
   
   % -- True Positive Ratio and the False Positive Ratio
   nb_CS=sum(powerband(:,1)); % Number of true CS
   nb_NCS=length(powerband(:,1))-nb_CS; % Number of true NCS
   nb_CS_ok=sum(res==3); % Number of good CS obtained
   nb_NCS_nok=sum(res==2); % Number of bad NCS obtained
   
   tpr=[nb_CS_ok/nb_CS, tpr]; % True Positive Ratio
   fpr=[nb_NCS_nok/nb_NCS, fpr]; % False Positive Ratio

end

%% FINAL THRESHOLD
dist=sqrt(fpr.^2+(1-tpr).^2); % Distance between top left corner and a point on the curve
[min_dist, argmin_dist]=min(dist); % Minimum distance 
final_threshold=thresholds_test(argmin_dist); % The better threshold

%% ACCURACY 
% -- Accuracy for all sections
final_label_obtained=(powerband(:,2)>final_threshold)*2; % 2=CS; 0=NCS
label_comparison=powerband(:,1)+final_label_obtained; % Comparison of true sections and the sections obtained (3=CS/CS; 1=CS/NCS; 2=NCS/CS; 0=NCS/NCS; )

% -- NCS accuracy
nb_section_NCS_OK=length(label_comparison(label_comparison==0)); % Number of NCS with a good label obtained 
nb_section_NCS=length(powerband(:,1)==0); % Total number of NCS
accuracy_NCS=nb_section_NCS_OK/nb_section_NCS *100; % Percentage of good labelling

% -- CS accuracy
nb_section_CS_OK=length(label_comparison(label_comparison==3)); % Number of CS with a good label obtained 
nb_section_CS=length(powerband(:,1)==1); % Total number of CS
accuracy_CS=nb_section_CS_OK/nb_section_CS *100; % Percentage of good labelling

% -- Total accuracy
accuracy=accuracy_NCS+accuracy_CS;

%% DISPLAY
figure, 
plot(fpr, tpr);
hold on 
plot(fpr, fpr, '--');
plot(fpr(argmin_dist),tpr(argmin_dist), '*');
hold off
title('ROC curve (with 500 thresholds)')
xlabel('False Positive Ratio')
ylabel('True Positive Ratio')
disp('OK thershold ROC');


% %% Autre méthode 
% cs=powerband(powerband(:,1)==1, 2);
% ncs=powerband(powerband(:,1)==0, 2);
% pred =[cs; ncs];

