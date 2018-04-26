% plotting subroutine
%
% julia dohner
% april 26, 2018

function myplot(residual10,delC10,yhat2,

figure
plot(residual10(:,1),residual10(:,2),delC10(:,1),yhat2)
xlabel('year')
ylabel('ppm CO2/year')
title('land uptake')
legend('Residual uptake','land uptake with T effects')
set(gca,'Xlim',[1850 2010])  
grid