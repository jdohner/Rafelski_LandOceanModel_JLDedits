% file BioBoxResponse.m
% 
% author Julia Dohner, borrowing from Lauren Rafelski
% 
% note Be sure to run defaults.m in outer folder before running this code
% 
% brief Calculates co2 uptake of fast land box model (driven by co2
% fertilization)
% 
% CO2 fertilization only
% code from bioboxtwo_sub10_annotate


function [C1dt,C2dt,delCdt,delC1,delC2] = BioBoxResponse(eps,Q1a,Q2a,ts,year,dpCO2a,T)

T0 = T(1,2); % initial temperature

dt = 1/ts; % timestep

% Define box sizes in ppm

Catm = 600/2.12; % around 283 ppm (preindustrial)

% fast biosphere box, from old model
% fast pool contains 110 PgC
C1 = 110/2.12; 
% slow biosphere box, changed to be same as 1 box
% slow pool contains 1477 PgC
C2 = 1477/2.12; 


% Rate constants - exchange in and out of fast box to atmosphere
 K1a = 1/2.5; 
 Ka1 = K1a*C1/Catm;

 % exchange constants in and out of slow box to atmosphere
K2a = 1/60; 
Ka2 = K2a*C2/Catm;

% set up arrays

% change in fast box
delC1(:,1) = year(:,1);
delC1(:,2) = zeros(size(year));
% change in slow box
delC2(:,1) = year(:,1);
delC2(:,2) = zeros(size(year));
% uptake by land?
delCdt(:,1) = year(:,1);

C1dt(:,1) = year(:,1);
C2dt(:,1) = year(:,1);

delC1(length(year)+1,1) = year(length(year),1)+dt;
delC2(length(year)+1,1) = year(length(year),1)+dt;

for i = 1:length(year)
    
    % fast box
    
    C1dt(i,2) = Ka1*(Catm + eps*dpCO2a(i,2)) - K1a*Q1a^((T(i,2)-T0)/10)*(C1 + delC1(i,2)); % temperature-dependent respiration
        
   % C1dt(i,2) = Ka1*(Catm + eps*dpCO2a(i,2))*(1 + Q1a*(T(i,2)-T0)) - K1a*(C1 + delC1(i,2)); % temperature-dependent photosynthesis
    
    % slow box
    
    C2dt(i,2) = Ka2*(Catm + eps*dpCO2a(i,2)) - K2a*Q2a^((T(i,2)-T0)/10)*(C2 + delC2(i,2)); % temperature-dependent respiration
        
   % C2dt(i,2) = Ka2*(Catm + eps*dpCO2a(i,2))*(1 + Q2a*(T(i,2)-T0)) - K2a*(C2 + delC2(i,2)); % temperature dependent photosynthesis  
    
    % box total change in concentrations
    
    delC1(i+1,2) = sum(C1dt(:,2))*dt;
    delC2(i+1,2) = sum(C2dt(:,2))*dt;
    
    % total flux into land
    
    delCdt(i,2) = C2dt(i,2) + C1dt(i,2);

    
end