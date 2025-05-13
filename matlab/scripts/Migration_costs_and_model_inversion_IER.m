global mu zstar dums MijH_data MijL_data U
global Nh_data Nl_data N_data Niih_data Niil_data 
global Lambda Gstar Hn Hr mH mL MijH MijL 
%-------------------------------------------------------------------------%
%                             Migration costs                             %
%-------------------------------------------------------------------------%
 cijH = max(0,1-(wH_data./wH_data').*((MijH_data./Niih_data).^mu));
 cijL = max(0,1-(wL_data./wL_data').*((MijL_data./Niil_data).^mu));


%-------------------------------------------------------------------------%
%                    Model inversion and exact mapping test               %
%-------------------------------------------------------------------------%
WijH = (wH_data'.*(1-cijH)).^(1/mu); 
WijL = (wL_data'.*(1-cijL)).^(1/mu); 

%% Implied demography and optimal location choice
SumWikH = wH_data.^(1/mu) + sum((1-dums).*WijH,1);
SumWikL = wL_data.^(1/mu) + sum((1-dums).*WijL,1);
mH = sum((1-dums).*WijH,1)./SumWikH;
mL = sum((1-dums).*WijL,1)./SumWikL; 
Niih = Nh_data.*((wH_data.^(1/mu))./SumWikH);
Niil = Nl_data.*((wL_data.^(1/mu))./SumWikL);


%% Skill premium
Lambda = SumWikH./SumWikL;

%% Exogenous supply of educational infrastructures 
Gstar = (Hn_data.^(1./(1+zstar)))./(1-1./Lambda);
Hn = (Gstar.*(1-1./Lambda)).^(1+zstar);
Nh = Hn.*N_data; Nl = (1-Hn).*N_data;

%% Endogeneized migration stock and resident human capital 
MijH = Nh.*(WijH./SumWikH);
MijL = Nl.*(WijL./SumWikL);
Lh = sum(MijH,2)'; Ll = sum(MijL,2)';
Ih = sum((1-dums).*MijH,2)'; Il = sum((1-dums).*MijL,2)';
Mh = sum((1-dums).*MijH,1); Ml = sum((1-dums).*MijL,1);
Hr = Lh./(Lh+Ll);



% Migration cost for validation exercise 
cijla= reshape(cijL,[U,1]);
cijha= reshape(cijH,[U,1]);
migH = migdta.Nijh; 
migL = migdta.Nijl;  
dumsa = migdta.dums;
o_iso = migdta.isoo;
d_iso = migdta.isod;
MIGCOST = table(o_iso,d_iso,dumsa,cijha,cijla,migH,migL);

% Export to CSV 
writetable(MIGCOST, 'output\MIGCOST.csv', 'Delimiter',',' ,'QuoteStrings', true)

