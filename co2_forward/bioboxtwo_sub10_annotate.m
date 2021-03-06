% calculate land uptake of co2 using co2 fertilization
% delCdt is total flux into land
%
% 12/18/07 - put in new turnover times for 110 Pg box
% 12/21/07 - change large box size to 1477 PgC, change corresponding atm
% flux
% 1/16/08 - change equations to allow for temperature dependent
% photosynthesis; still call it Q so it's easy to run

function [C1dt,C2dt,delCdt,delC1,delC2] = bioboxtwo_sub10_annotate(eps,Q1a,Q2a,ts,year,dpCO2a,T)

T0 = T(1,2);

dt = 1/ts;

% Define box sizes in ppm

Catm = 600/2.12; % around 283 ppm (preindustrial)

C1 = 110/2.12; % fast biosphere box, from old model
C2 = 1477/2.12; % slow biosphere box, changed to be same as 1 box


% Rate constants
 K1a = 1/2.5; 
 Ka1 = K1a*C1/Catm;

K2a = 1/60; % slow box to atmosphere
Ka2 = K2a*C2/Catm;

% set up arrays

delC1 = [year(:,1) zeros(size(year))];
delC2 = [year(:,1) zeros(size(year))];
% delC1(:,1) = year(:,1);
% delC1(:,2) = zeros(size(year));
% delC2(:,1) = year(:,1);
% delC2(:,2) = zeros(size(year));
delCdt(:,1) = year(:,1);
C1dt(:,1) = year(:,1);
C2dt(:,1) = year(:,1);


delC1(length(year)+1,1) = year(length(year),1)+dt;
delC2(length(year)+1,1) = year(length(year),1)+dt;

% this loop is the same as in oceanpulseresponse- need to be the same as in
% ocean pulse response in higher level code, then call next what's below
for ii = 1:length(year)
    
    % fast box (equation (3) in Rafelski 2009
    % calculate land uptake by fast box
    C1dt(ii,2) = Ka1*(Catm + eps*dpCO2a(ii,2)) - K1a*Q1a^((T(ii,2)-T0)/10)*(C1 + delC1(ii,2)); % temperature-dependent respiration
        
   % C1dt(i,2) = Ka1*(Catm + eps*dpCO2a(i,2))*(1 + Q1a*(T(i,2)-T0)) - K1a*(C1 + delC1(i,2)); % temperature-dependent photosynthesis
    
    % slow box (equation (3) in Rafelski 2009
    % calculate land uptake by slow box
    C2dt(ii,2) = Ka2*(Catm + eps*dpCO2a(ii,2)) - K2a*Q2a^((T(ii,2)-T0)/10)*(C2 + delC2(ii,2)); % temperature-dependent respiration
        
   % C2dt(i,2) = Ka2*(Catm + eps*dpCO2a(i,2))*(1 + Q2a*(T(i,2)-T0)) - K2a*(C2 + delC2(i,2)); % temperature dependent photosynthesis  
    
    % box total change in concentrations (of fast box, slow box, respectively)
    delC1(ii+1,2) = sum(C1dt(:,2))*dt;
    delC2(ii+1,2) = sum(C2dt(:,2))*dt;
    
    % total flux into land
    
    delCdt(ii,2) = C2dt(ii,2) + C1dt(ii,2);
    
end



