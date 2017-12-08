% savedPlot = plot(residualLandUptake(:,1),residualLandUptake(:,2),'g',ff1(:,1),ff1(:,2),'k',dtdelpCO2a_obs(:,1),dtdelpCO2a_obs(:,2),'r',fas(:,1),-Aoc*fas(:,2),'b');
% axis(savedPlot,[1800 2010 -10 10])
% legend(savedPlot,'residualLandUptake','fossil fuel','atmosphere','ocean','Location','SouthWest')
% title(savedPlot,'residualLandUptake plus components JLD ')
% xlabel(savedPlot,'Year ')
% ylabel(savedPlot,'ppm/year')
% 
% savefig(savedPlot,'JLD_plot.fig');
% 
% h1 = openfig('LR_plot.fig','reuse'); % open figure
% ax1 = gca; % get handle to axes of figure
% h2 = openfig('JLD_plot.fig','reuse');
% ax2 = gca;
% h3 = figure;
% s1 = subplot(2,1,1); %create and get handle to the subplot axes
% s2 = subplot(2,1,2);
% fig1 = get(ax1,'children'); %get handle to all the children in the figure
% fig2 = get(ax2,'children');
% copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
% copyobj(fig2,s2);
% 
% % axis(s1,[1800 2010 -10 10])
% % legend(s1,'fossil fuel','atmosphere','ocean','land','Location','SouthWest')
% % title(s1,'Sources and sinks from Joos response function LR')
% % xlabel(s1,'Year ')
% % ylabel(s1,'ppm/year')
% 
% % axis(s2,[1800 2010 -10 10])
% % legend(s2,'residualLandUptake','fossil fuel','atmosphere','ocean','Location','SouthWest')
% % title(s2,'residualLandUptake plus components JLD ')
% % xlabel(s2,'Year ')
% % ylabel(s2,'ppm/year')


% %% suplot option #2:
% 
% % Load saved figures
% c=hgload('LR_plot.fig');
% ax1 = gca;
% k = plot(residualLandUptake(:,1),residualLandUptake(:,2),'g',ff1(:,1),ff1(:,2),'k',dtdelpCO2a_obs(:,1),dtdelpCO2a_obs(:,2),'r',fas(:,1),-Aoc*fas(:,2),'b');
% ax2 = gca;
% %k=hgload('MySecondFigure.fig');
% % Prepare subplots
% figure
% h(1)=subplot(2,1,1);
% h(2)=subplot(2,1,2);
% % Paste figures on the subplots
% copyobj(allchild(get(c,'CurrentAxes')),h(1));
% copyobj(allchild(get(k,'CurrentAxes')),h(2));
% % Add legends
% l(1)=legend(h(1),'LegendForFirstFigure')
% l(2)=legend(h(2),'LegendForSecondFigure')