global sigma psi eta U
global dums L_data Hn Hr MijH MijL 
%-------------------------------------------------------------------------%
%                         Initial Equilibrium                             %
%-------------------------------------------------------------------------%

 %          Development accounting 
WR = wH_data./wL_data;
SRL = Hr./(1-Hr);
thetH = WR./(WR+SRL.^(-1/sigma));
thetL = 1-thetH;
LiT = (thetH.*Hr.^psi + thetL.*(1-Hr).^psi).^(1/psi);
w = Hr.*wH_data+(1-Hr).*wL_data;
TFP = w./LiT;

  % wage replication 
  wH = TFP.*thetH.*(LiT./Hr).^(1/sigma);
  wL = TFP.*thetL.*(LiT./(1-Hr)).^(1/sigma); 

  % Disposable income (income adjusted to fiscal burden)
  mij = (MijH+MijL)./(Nh+Nl);
  mi = sum((1-dums).*mij,1); % Outmigration rate 
  tau1 = GovCons_data.*L_data.^(-eta) + EduExp_data./(1-mi);
  ywd  =  w.*(1-tau1); 

 %% Remittances 
% Pseudo endogeneization of remittances (just to let them endogeneously
% adjust with the size of migration. So no behavioural change here) 

    % Relative size of corridors 
    Mij_data = (1-dums).*MijH_data+(1-dums).*MijL_data;
    Mi_data = sum(Mij_data,1);
    wgt = Mij_data./Mi_data; 
    Rij_data = RemiT_data.*wgt; % Spread remittances across corridors
    sum(Rij_data(:)) % $472.24Bn 
    sum(RemiT_data(:)) % $472.24Bn 
    Rij_pm = Rij_data./Mij_data; Rij_pm(isnan(Rij_pm))=0; % remittances per migrant across corridors 
    Rij = Rij_pm.*Mij_data;
    r_data = sum(Rij,1)./L_data; % relative size of remittances gain in percentage of per worker disposable income

%% Income per natural (Clemens and Pritchett)
  wijn = (wL'.*(1-tau1').*MijL + wH'.*(1-tau1').*MijH);
  ypnd = sum(wijn,1)./N_data;  % Income per natural in absence of migration 

  wijnm = (wL.*(1-tau1).*Nl_data + wH.*(1-tau1).*Nh_data);
  yprd  = sum(wijnm,1)./N_data;  % Income per natural in absence of migration 
  dwn   = ypnd./yprd; % relative deviation in income per natural due to migration 
  

 %% Income inequality (Theil index) 
 wHd = wH.*(1-tau1);
 wLd = wL.*(1-tau1);
 wd  = Hr.*wHd + (1-Hr).*wLd;
 wdbar = sum(wd.*L_data)./sum(L_data); % World average disposable income 
 Theil = (sum(wHd.*Lh_data.*log(wHd./wdbar)) + sum(wLd.*Ll_data.*log(wLd./wdbar)))./(wdbar.*sum(L_data)); % Theil index 
 TheilA = (sum(wd.*L_data.*log(wd./wdbar)))./(wdbar.*sum(L_data)); % Theil index across countries 
 InnerT = (((wHd.*Lh_data).*log(wHd./wd)) + ((wLd.*Ll_data).*log(wLd./wd)))./(wd.*L_data); % inner component of the across 
 TheilW = sum((wd.*L_data.*InnerT))./(wdbar.*sum(L_data)); % Theil index across countries 

