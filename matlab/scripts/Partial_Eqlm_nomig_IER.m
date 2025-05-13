global mu psi kappa epsilon eta lambda1 rho mbar o d zstar dums  
global N_data L_data Gstar Hr 
%-------------------------------------------------------------------------%
%                       Partial Equilibrium (nomig)                       %
%-------------------------------------------------------------------------%

cijHNew = (1-dums).*ones(o,d); % Shock  
cijLNew = (1-dums).*ones(o,d); % Shock  

WijHNew = (wH'.*(1-cijHNew)).^(1/mu); % Wages do not respond endogeneously here 
WijLNew = (wL'.*(1-cijLNew)).^(1/mu); 

%% Implied demography
SumWijHNew =  sum(WijHNew,1);
SumWijLNew =  sum(WijLNew,1);

SumWikHNew = wH.^(1/mu) + sum((1-dums).*WijHNew,1);
SumWikLNew = wL.^(1/mu) + sum((1-dums).*WijLNew,1);

mHNew = sum((1-dums).*WijHNew,1)./SumWikHNew;
mLNew = sum((1-dums).*WijLNew,1)./SumWikLNew; 
LambdaNew = SumWikHNew./SumWikLNew;
HnNew = (Gstar.*(1-1./LambdaNew)).^(1+zstar);
NhNew = HnNew.*N_data; NlNew = (1-HnNew).*N_data;
MijHNew = NhNew.*(WijHNew./SumWikHNew);
MijLNew = NlNew.*(WijLNew./SumWikLNew);
LhNew = sum(MijHNew,2)'; LlNew = sum(MijLNew,2)';
IhNew = sum((1-dums).*MijHNew,2)'; Il = sum((1-dums).*MijLNew,2)';
MhNew = sum((1-dums).*MijHNew,1); MlNew = sum((1-dums).*MijLNew,1);
HrNew = LhNew./(LhNew+LlNew);
NiihNew = NhNew.*((wH.^(1/mu))./SumWikHNew);
NiilNew = NlNew.*((wL.^(1/mu))./SumWikLNew);


 %-------------------------------------------------------------------------
 %          Development accounting (nomig)
 %-------------------------------------------------------------------------
LiTNew = (thetH.*HrNew.^psi + thetL.*(1-HrNew).^psi).^(1/psi);
ywdNew =TFP.*LiTNew.*(1-tau1);
NIR1_NM  = ywd./ywdNew -1; % Human capital channel

SRN = HrNew./(1-HrNew);
R1 = SRN./SRL;
R_Thet = thetH./thetL;
thetHNew =(R_Thet.*R1.^kappa)./(1+(R_Thet.*R1.^kappa));
thetLNew =  1-thetHNew;
LiTNew = (thetHNew.*HrNew.^psi + thetLNew.*(1-HrNew).^psi).^(1/psi);
TFPNew = TFP.*R1.^epsilon;
ywdNew =TFPNew.*LiTNew.*(1-tau1);
NIR2_NM = ywd./ywdNew -1;

denomi = (thetHNew.*Hr.^psi +thetLNew.*(1-Hr).^psi).^(1/psi); % technological externality only
numera = (thetH.*Hr.^psi + thetL.*((1-Hr).^psi)).^(1/psi);
TEC_NM  = (TFP./TFPNew - 1) +(numera./denomi -1); 

 TFPNew  = TFPNew.*((mbar./(mbar+ mi)).^rho); % Diaspora externalities
 ywdNew  = TFPNew.*LiTNew.*(1-tau1);
 NIR3_NM = ywd./ywdNew -1;

 TFPNew = TFP.*(mbar./(mbar+ mi)).^rho; % Diaspora externality only
 DIAS_NM  = TFP./TFPNew - 1; 

  LNew = LhNew+LlNew;
  mijNew = (MijHNew+MijLNew)./(NhNew+NlNew);
  miNew = sum((1-dums).*mijNew,1); % Outmigration rate New 
  tau1_New = GovCons_data.*LNew.^(-eta) + EduExp_data./(1-miNew); % Account for fiscal externalities (tax and size)
  ywdNew  = TFPNew.*LiTNew.*(1-tau1_New);
  NIR4_NM = NIR3_NM + (tau1_New - tau1)./(1-tau1_New);   

  FIS_NM  = (tau1_New - tau1)./(1-tau1_New); % fiscal externality only

  NQ =LiT ; % Account for market size effects
  NQNew  = NQ./(1-mi);     LossNQ =(NQ - NQNew)./NQNew;
  NIR5_NM  = NIR4_NM + 0.7*(LossNQ./(lambda1 - 1));
  MKT_NM  = 0.7*(LossNQ./(lambda1 - 1)); % market externality only
    
  ri = RemiT_data./L_data; %   Account for remittances
  ywd1   = ywd + ri;
  REM_NM  =ri./ywdNew;  % remittances externality only
  NIR6_NM  = NIR5_NM + REM_NM; %ywd1./ywdNew-1 + 0.7*(LossNQ./(lambda1 - 1));
  ywd_nomig= ywdNew.*(1-MKT_NM);
 

% Income per natural
    wijNew = dums.*(wL.*(1-tau1_New).*NiilNew + wH.*(1-tau1_New).*NiihNew) + (1-dums).*(wL'.*(1-tau1').*MijLNew + wH'.*(1-tau1').*MijHNew);
    ypndNew = sum(wijNew,1)./N_data; 
    dwnNew   = ypnd./ypndNew-1; % relative deviation in income per natural due to migration 
    NIR7_NM = NIR5_NM+dwnNew; % remittances are factored out as already part of migrants' income abroad
    NIR7_NM = min(.9,NIR7_NM);


%% World averages (NIR6, NIR7)
LNM = LhNew+LlNew;
NIR6_NMworld = sum(NIR6_NM.*LNM)./sum(LNM); % World average net disposable income response 
NIR7_NMworld = sum(NIR7_NM.*LNM)./sum(LNM); % World average net disposable income response 

%% Low income countries averages  (NIR6, NIR7)
NIR6_NMLIC = sum(LIC.*NIR6_NM.*LNM)./sum(LIC.*LNM); 
NIR7_NMLIC = sum(LIC.*NIR7_NM.*LNM)./sum(LIC.*LNM); 

%% Lower middle income countries averages  (NIR6, NIR7)
NIR6_NMLMIC = sum(LMIC.*NIR6_NM.*LNM)./sum(LMIC.*LNM); 
NIR7_NMLMIC = sum(LMIC.*NIR7_NM.*LNM)./sum(LMIC.*LNM); 

%% Upper middle income countries averages  (NIR6, NIR7)
NIR6_NMUMIC = sum(UMIC.*NIR6_NM.*LNM)./sum(UMIC.*LNM); 
NIR7_NMUMIC = sum(UMIC.*NIR7_NM.*LNM)./sum(UMIC.*LNM); 

%% High income countries averages  (NIR6, NIR7)
NIR6_NMHIC = sum(HIC.*NIR6_NM.*LNM)./sum(HIC.*LNM); 
NIR7_NMHIC = sum(HIC.*NIR7_NM.*LNM)./sum(HIC.*LNM); 

%--------------------------------------------------------------------------
%   CHANNELS
%--------------------------------------------------------------------------
LNM = LNew; LhNM= LhNew; LlNM= LlNew; HrNM = HrNew; HnNM = HnNew; LambdaNM = LambdaNew;
dLamb = Lambda./LambdaNM;
c_HumCap = NIR1_NM;
c_TecExt = NIR2_NM - NIR1_NM;
c_DiaExt = NIR3_NM - NIR2_NM;
c_FisExt = NIR4_NM - NIR3_NM;
c_MktExt = NIR5_NM - NIR4_NM;
c_RemEff = NIR6_NM - NIR5_NM;
c_NatEff = NIR7_NM - NIR6_NM;

cc_HumCap = NIR1_NM;
cc_TecExt = TEC_NM;
cc_DiaExt = DIAS_NM;
cc_FisExt = FIS_NM;
cc_MktExt = MKT_NM;
cc_RemEff = REM_NM;
cc_NatEff = NIR7_NM - NIR6_NM;

c_channels =[c_HumCap',c_TecExt',c_DiaExt',c_FisExt',c_MktExt',c_RemEff',c_NatEff',cc_HumCap',cc_TecExt',cc_DiaExt',cc_FisExt',cc_MktExt',cc_RemEff',cc_NatEff',NIR6_NM',NIR7_NM',ywd_nomig',ywd1'];
cols = ["HumCap","TecExt","DiaExt","FisExt","MktExt","RemEff","NatEff","HumCap2","TecExt2","DiaExt2","FisExt2","MktExt2","RemEff2","NatEff2","NIR6_NMb","NIR7_NMb","ynomig_b","ywd"];
cc_channels =[cols;c_channels];
