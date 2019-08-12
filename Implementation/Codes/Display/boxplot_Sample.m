clear
close all
clc
F1=randn(40,4);
F2=randn(40,4).*2;
%F3=randn(40,4).*2;

position1=1:size(F1,2);
position2=position1+0.15;

figure;
boxplot(F1,'Color', 'b','positions', position1,'width',0.12)
hold on
boxplot(F2,'Color', 'r','positions', position2,'width',0.12)
