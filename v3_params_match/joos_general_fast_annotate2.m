% Pulse response model using response function "r"
% 7/15/09: add "dpCO2s" to output

function [fas,dpCO2s,T]= joos_general_fast_annotate2(yearOcean,dpCO2a,c,h,kg,Tconst,Aoc,r,dt,varSST)

dpCO2s = zeros(length(dpCO2a),2); % dissolved CO2
dpCO2s(:,1) = dpCO2a(:,1);
fas = zeros(length(yearOcean),2);
integral = zeros(length(yearOcean),length(yearOcean));
delDIC = zeros(length(yearOcean),2);

if varSST == 1
    load newsst_2011.mat;
    hadcrut = csvread('HadCRUT4.csv');
    % interpolate to monthly
    sstmo(:,1) = yearOcean;
    sstmo(:,2) = interp1(hadcrut(1:2:end,1),hadcrut(1:2:end,14),yearOcean);
    i = find(sstmo(:,1) >= 1850,1);
    sstmo(1:i-1,2) = 0; % fill in values before 1850
    %sstmo(:,2) = sstmo(:,2) + meanSST;
    T = sstmo;
    
    figure
    plot(T(:,1),T(:,2));
    xlabel('year');
    ylabel('degree C');
    title('sst anomaly record');
end

for m = 1:(length(yearOcean)-1)
    %Calculate flux 
   fas(m,1) = yearOcean(m);
   fas(m,2) = (kg/Aoc)*(dpCO2a(m,2) - dpCO2s(m,2)); % air-sea flux of CO2
   
    w = conv(fas(1:m,2),r(1:m,2)); % convolve the air-sea flux and the pulse response function, as in Joos 1996

    % Calculate delDIC 
    delDIC(m+1,1) = yearOcean(m+1);
   
    delDIC(m+1,2) = (c/h)*w(m)*dt; % change in DIC

    %Calculate dpCO2s from DIC - from Joos 1996
    if varSST == 0 % fixed sst
        dpCO2s(m+1,2) = (1.5568 - (1.3993E-2)*Tconst)*delDIC(m+1,2) + (7.4706-0.20207*Tconst)*10^(-3)*...
        (delDIC(m+1,2))^2 - (1.2748-0.12015*Tconst)*10^(-5)*(delDIC(m+1,2))^3 + (2.4491-0.12639*Tconst)...
        *10^(-7)*(delDIC(m+1,2))^4 - (1.5468-0.15326*Tconst)*10^(-10)*(delDIC(m+1,2))^5;
    
    else % variable sst
        dpCO2s(m+1,2) = (1.5568 - (1.3993E-2)*T(m+1,2))*delDIC(m+1,2) + (7.4706-0.20207*T(m+1,2))*10^(-3)*...
        (delDIC(m+1,2))^2 - (1.2748-0.12015*T(m+1,2))*10^(-5)*(delDIC(m+1,2))^3 + (2.4491-0.12639*T(m+1,2))...
        *10^(-7)*(delDIC(m+1,2))^4 - (1.5468-0.15326*T(m+1,2))*10^(-10)*(delDIC(m+1,2))^5;
    
        %dpCO2s(m+1,2) = dpCO2s(m+1,2)*exp(0.0423*T(m+1,2));
    end
end

% calculate the flux for the last time point
fas(length(yearOcean),1) = yearOcean(length(yearOcean));
fas(length(yearOcean),2) = (kg/Aoc)*(dpCO2a(length(yearOcean),2) - dpCO2s(length(yearOcean),2));

zero1 = find(delDIC(:,1) == 0);
delDIC(zero1,2) = NaN;
delDIC(zero1,1) = NaN;

zero2 = find(dpCO2s(:,2) == 0);
dpCO2s(zero2,2)  = NaN;

if varSST == 0
    T = Tconst;
end
