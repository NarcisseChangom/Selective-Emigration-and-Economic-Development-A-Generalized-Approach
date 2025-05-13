global mu psi kappa epsilon eta lambda1 rho mbar zstar 
global N_data L_data Niih Niil nijH nijL nijbar
global Gstar Hr Hn MijH MijL 

%-------------------------------------------------------------------------%
%                       Partial Equilibrium (noselec)                       %
%-------------------------------------------------------------------------%

%% movers to stayers ratios 
nijH = MijH./Niih;
nijL = MijL./Niil;
nijbar = Hn.*nijH + (1-Hn).*nijL;

%% migration costs compatible with neutral selection 
cijH_NS = 1 - (wH./wH').*nijbar.^mu;
cijL_NS = 1 - (wL./wL').*nijbar.^mu;

%% test 
nijH_hat = ((wH'./wH).*(1-cijH_NS)).^(1/mu); % It lowers the prob of HS and enlarges the prob of LS
nijL_hat = ((wL'./wL).*(1-cijL_NS)).^(1/mu);
scatter((1-dums).*cijH,(1-dums).*cijH_NS)
scatter((1-dums).*cijL,(1-dums).*cijL_NS)
scatter(nijH,nijH_hat)
scatter(nijL,nijL_hat)

%% Implied counterfactual demography and spatial re-allocation of workers

WijHNew = (wH'.*(1-cijH_NS)).^(1/mu); % Wages do not respond endogeneously here 
WijLNew = (wL'.*(1-cijL_NS)).^(1/mu); 
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
 %          Development accounting (noselec)
 %-------------------------------------------------------------------------
LiTNew = (thetH.*HrNew.^psi + thetL.*(1-HrNew).^psi).^(1/psi);
ywdNew =TFP.*LiTNew.*(1-tau1);
NIR1_NS  = ywd./ywdNew -1; % Human capital channel

SRN = HrNew./(1-HrNew);
R1 = SRN./SRL;
R_Thet = thetH./thetL;
thetHNew =(R_Thet.*R1.^kappa)./(1+(R_Thet.*R1.^kappa));
thetLNew =  1-thetHNew;
LiTNew = (thetHNew.*HrNew.^psi + thetLNew.*(1-HrNew).^psi).^(1/psi);
TFPNew = TFP.*R1.^epsilon;
ywdNew =TFPNew.*LiTNew.*(1-tau1);
NIR2_NS = ywd./ywdNew -1;

denomi = (thetHNew.*Hr.^psi +thetLNew.*(1-Hr).^psi).^(1/psi); % technological externality only
numera = (thetH.*Hr.^psi + thetL.*((1-Hr).^psi)).^(1/psi);
TEC_NS  = (TFP./TFPNew - 1) +(numera./denomi -1); 

 TFPNew = TFPNew.*((mbar./(mbar+ mi)).^rho); % Diaspora externalities
 ywdNew  = TFPNew.*LiTNew.*(1-tau1);
 NIR3_NS= ywd./ywdNew -1;

 TFPNew = TFP.*(mbar./(mbar+ mi)).^rho; % Diaspora externality only
 DIAS_NS  = TFP./TFPNew - 1; 

  LNew = LhNew+LlNew;
  mijNew = (MijHNew+MijLNew)./(NhNew+NlNew);
  miNew = sum((1-dums).*mijNew,1); % Outmigration rate New 
  tau1_New = GovCons_data.*LNew.^(-eta) + EduExp_data./(1-miNew); % Account for fiscal externalities (tax and size)
  ywdNew  = TFPNew.*LiTNew.*(1-tau1_New); 
  NIR4_NS    = ywd./ywdNew - 1; 
  FIS_NS  = (tau1_New - tau1)./(1-tau1_New); % fiscal externality only

  NQ =LiT ; % Account for market size effects
  NQNew  = NQ.*((1-miNew)./(1-mi));
  LossNQ =(NQ - NQNew)./NQNew;
  NIR5_NS  = NIR4_NS + 0.7*(LossNQ./(lambda1 - 1));
  MKT_NS  = 0.7*(LossNQ./(lambda1 - 1)); % market externality only

  ri_NS = RemiT_data./(Nii_data);
  REM_NS  = ri_NS./ywdNew;  % remittances externality only
  ywd2 = ywd+r_data;
  ywdNew = ywdNew+ri_NS;

  NIR6_NS  = ywd2./ywdNew-1 + MKT_NS; % 0.7*(LossNQ./(lambda1 - 1));
  ywd_noselec= ywdNew.*(1-MKT_NS);


% Income per natural
    %wijNew = dums.*(wL.*(1-tau1_New).*NiilNew + wH.*(1-tau1_New).*NiihNew) + (1-dums).*(wL'.*(1-tau1').*MijLNew + wH'.*(1-tau1').*MijHNew);
    wijNew = dums.*(wL.*(1-tau1_New).*NlNew + wH.*(1-tau1_New).*NhNew); %+ (1-dums).*(wL'.*(1-tau1').*MijLNew + wH'.*(1-tau1').*MijHNew);
    ypndNew = sum(wijNew,1)./N_data; 
    dwnNew   = ypnd./ypndNew-1; % relative deviation in income per natural due to migration 
    NIR7_NS = NIR5_NS+dwnNew; % remittances are factored out as already part of migrants' income abroad
    NIR7_NS = min(.9,NIR7_NS);


%% World averages (NIR6, NIR7)
LNS = LhNew+LlNew;
NIR6_NSworld = sum(NIR6_NS.*LNS)./sum(LNS); % World average net disposable income response 
NIR7_NSworld = sum(NIR7_NS.*LNS)./sum(LNS); % World average net disposable income response 

%% Low income countries averages  (NIR6, NIR7)
NIR6_NSLIC = sum(LIC.*NIR6_NS.*LNS)./sum(LIC.*LNS); 
NIR7_NSLIC = sum(LIC.*NIR7_NS.*LNS)./sum(LIC.*LNS); 

%% Lower middle income countries averages  (NIR6, NIR7)
NIR6_NSLMIC = sum(LMIC.*NIR6_NS.*LNS)./sum(LMIC.*LNS); 
NIR7_NSLMIC = sum(LMIC.*NIR7_NS.*LNS)./sum(LMIC.*LNS); 

%% Upper middle income countries averages  (NIR6, NIR7)
NIR6_NSUMIC = sum(UMIC.*NIR6_NS.*LNS)./sum(UMIC.*LNS); 
NIR7_NSUMIC = sum(UMIC.*NIR7_NS.*LNS)./sum(UMIC.*LNS); 

%% High income countries averages  (NIR6, NIR7)
NIR6_NSHIC = sum(HIC.*NIR6_NS.*LNS)./sum(HIC.*LNS); 
NIR7_NSHIC = sum(HIC.*NIR7_NS.*LNS)./sum(HIC.*LNS);

%--------------------------------------------------------------------------
%   CHANNELS
%--------------------------------------------------------------------------
LNS = LNew; LhNS= LhNew; LlNS= LlNew; HrNS = HrNew; HnNS = HnNew; LambdaNS = LambdaNew;
c_HumCap = NIR1_NS;
c_TecExt = NIR2_NS - NIR1_NS;
c_DiaExt = NIR3_NS - NIR2_NS;
c_FisExt = NIR4_NS - NIR3_NS;
c_MktExt = NIR5_NS - NIR4_NS;
c_RemEff = NIR6_NS - NIR5_NS;
c_NatEff = NIR7_NS - NIR6_NS;

cc_HumCap = NIR1_NS;
cc_TecExt = TEC_NS;
cc_DiaExt = DIAS_NS;
cc_FisExt = FIS_NS;
cc_MktExt = MKT_NS;
cc_RemEff = REM_NS;
cc_NatEff = NIR7_NS - NIR6_NS;

c_channels =[c_HumCap',c_TecExt',c_DiaExt',c_FisExt',c_MktExt',c_RemEff',c_NatEff',cc_HumCap',cc_TecExt',cc_DiaExt',cc_FisExt',cc_MktExt',cc_RemEff',cc_NatEff',NIR6_NS',NIR7_NS',ywd_noselec',ywd1'];
cols = ["HumCap","TecExt","DiaExt","FisExt","MktExt","RemEff","NatEff","HumCap2","TecExt2","DiaExt2","FisExt2","MktExt2","RemEff2","NatEff2","NIR6_NSb","NIR7_NSb","ynoselec_b","ywd"];
cc_channels_NS =[cols;c_channels];
