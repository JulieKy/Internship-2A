% PLOT: Test the plots

x = 0:pi/100:2*pi;
y = sin(x);

figure,
hax=axes;

band_end=length(x);
x_interval=floor(length(x)/4)

for n_band = 1 : floor(length(x)/x_interval)
    band_start=band_end-(floor(x_interval));
    h_axe=get(hax,'YLim');
    
    %% Rectangle
    x_rec=x(band_start);
    y_rec=h_axe(1);
    w_rec=x_interval;
    h_rec=h_axe(2)-h_axe(1);
    rectangle('Position',[x(band_start) y_rec 2 h_rec],'FaceColor',[0.5 .5 .5],'EdgeColor','b',...
        'LineWidth',1); 
    hold on
    
    %% Lines
    line([x(band_start) x(band_start)],h_axe, 'Color',[0 0 0]); % Vertical lines differentiating the frequency bands
    line([x(band_start) x(band_end)], [mean(y(band_start:band_end)) mean(y(band_start:band_end))], 'LineWidth',2, 'Color',[1 0 0]); % Horizontal lines (mean)
    
    band_end=band_start;
    
end
 
%% Signal 
plot(x,y, 'LineWidth',4, 'Color',[0.4 1 0.4]); 

hold off
