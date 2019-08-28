% PLOT: Test the plots

% x = 0:pi/100:2*pi;
% y = sin(x);
% 
% % figure,
% % hax=axes;
% % 
% % band_end=length(x);
% % x_interval=floor(length(x)/4)
% % 
% % % for n_band = 1 : floor(length(x)/x_interval)
% %     band_start=band_end-(floor(x_interval));
% %     h_axe=get(hax,'YLim');
% %     
% %     %% Rectangle
% %     x_rec=x(band_start);
% %     y_rec=h_axe(1);
% %     w_rec=x_interval;
% %     h_rec=h_axe(2)-h_axe(1);
% %     rectangle('Position',[x(band_start) y_rec 2 h_rec],'FaceColor',[0.5 .5 .5],'EdgeColor','b',...
% %         'LineWidth',1); 
% %     hold on
% %     
% %     %% Lines
% %     line([x(band_start) x(band_start)],h_axe, 'Color',[0 0 0]); % Vertical lines differentiating the frequency bands
% %     line([x(band_start) x(band_end)], [mean(y(band_start:band_end)) mean(y(band_start:band_end))], 'LineWidth',2, 'Color',[1 0 0]); % Horizontal lines (mean)
% %     
% %     band_end=band_start;
% %     
% % end
% %  
% % %% Signal 
% % plot(x,y, 'LineWidth',4, 'Color',[0.4 1 0.4]); 
% % 
% % hold off
% 
% 
% %% Box plot 
% % figure, 
% y1=[2 3 4 3 4 5]';
% y2=[4 5 6 7 6 8]'; 
% y3=[2 9 1 3 4 5]';
% y4=[5 6 6 6 6 7]'; 
% y6=[5 6 6 6 6 7]'; 
% y5=[5 6 6 6 6 7]'; 
% y7=[5 6 6 6 6 7]'; 
% 
% % boxplot([y1,y2],'Labels',{'mu = 5','mu = 6'})
% % title('Boxplot test')
% % xlabel('Country of Origin')
% % ylabel('Miles per Gallon (MPG)')
% 
% % figure
% % boxplot([y1,y2, y4, y3, y5, y6, y7],'Labels',{'mu = 5', 'mu = 5','mu = 6', 'mu = 7','mu = 8', 'mu = 6', 'mu = 7','mu = 8'}); 
% % dim = [.2 .5 .3 .3];
% % str = 'Straight Line Plot from 1 to 10';
% % annotation('textbox',dim,'String',str,'FitBoxToText','on');
% % title('Compare Random Data from Different Distributions')
% % dim = [0 0.4 0.2 0.6];
% % annotation('rectangle',dim,'Color','red')
% % % load carsmall
% % figure,
% % boxplot(MPG,Origin)
% % title('Miles per Gallon by Vehicle Origin')
% % xlabel('Country of Origin')
% % ylabel('Miles per Gallon (MPG)')
% 
% 
% figure, 
% title('truc'); 
% subplot(2,1,1); plot(y2,y1) title('truc'); 
% subplot(2,1,2); plot(y1,y2)
% 
% 
% 
 length_time=randn(1,50)*100; min_ok=10;

figure, 
hax=axes;
h_axe=get(hax,'XLim');
length_short=length_time(length_time<min_ok);
length_ok=length_time(length_time>=min_ok);
plot(length_short, '*', 'Color',[0.8 0 0]); hold on
plot(length_ok, '*', 'Color',[0 0.6 0]);
plot(1:length(length_time), length_time, '*'); hold on
%plot(1:length(length_time),ones(1,length(length_time))*10 , '--', 'Color', 'r');
line([1,length(length_time)],[10,10], 'Color',[0.8 0 0], 'LineWidth', 2);
plot(10, 11, '-o', 'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6])
text(10,11,'\leftarrow sin(\pi)')
hold off
title('Duration of Samples after Cry Removal')
xlabel('Samples')
ylabel('Duration [s]')
