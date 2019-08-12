% TEST K-FOLD CROSS VALIDATION

% Load the data set
load fisheriris

% Create indices for the 10-fold cross-validation
k=6;
indices = crossvalind('Kfold',species,k); % Number of their fold

% Initialize an object to measure the performance of the classifier
cp = classperf(species);

% Perform the classification using the measurement data and report the error rate
for i = 1:10
    test = (indices == i);  % One fold for Validation
    train = ~test; % k-1 folds for training
    class = classify(meas(test,:),meas(train,:),species(train,:));
    classperf(cp,class,test);
end
EQ=cp.ErrorRate % The ratio of the number of incorrectly classified samples divided by the total number of classified samples
