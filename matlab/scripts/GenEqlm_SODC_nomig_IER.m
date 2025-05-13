

%-------------------------------------------------------------------------%
%                GENERAL EQUILIBRIUM: SODC (ENDO. WAGES  ORIG.)           %
%-------------------------------------------------------------------------%
% No migration 
% SODC = Small Openned Developing Country 

global sigma psi mu kappa epsilon eta lambda1 rho mbar zstar thetHPE thetLPE
global N_data L_data Gstar HrPE HnPE LambdaPE LiTPE  o d dums TFPPE tauPE 

%------------------------
% A: Equilibrium replication  
%------------------------ 
cijHPE = (1-dums).*ones(o,d); cijLPE =  (1-dums).*ones(o,d); % cijHPE = cijH; cijLPE = cijL;
%cijHPE = cijH; cijLPE = cijL;
LiTPE = LiT; thetHPE = thetH; thetLPE = 1-thetHPE; TFPPE = TFP; tauPE = tau1;
LhPEtemp = Lh; LlPEtemp = Ll;


corr=0.9; maxerror =1.0e-5; maxit = 5000; error=2; count = 0; 
 while (error > maxerror) && (count < maxit)
      count = count + 1;
       error=0; errora=0;  errorb=0;

    % Skill specific endogeneous wages 
    LPEtemp = LhPEtemp+LlPEtemp;
    HrPE = LhPEtemp./LPEtemp;
    wHPE = thetHPE.*TFPPE.*(LiTPE./HrPE).^(1/sigma);
    wLPE = thetLPE.*TFPPE.*(LiTPE./(1-HrPE)).^(1/sigma);

    % Corridor specific utility and optimal re-allocation choice 
    WijHPE = (wH'.*(1-cijHPE)).^(1/mu);  
    WijLPE = (wL'.*(1-cijLPE)).^(1/mu); 

    SumWikHPE = wHPE.^(1/mu) + sum((1-dums).*WijHPE,1); % endogeneous wages at origin and exogenous at destination 
    SumWikLPE = wLPE.^(1/mu) + sum((1-dums).*WijLPE,1);

    mHPE = sum((1-dums).*WijHPE,1)./SumWikHPE;
    mLPE = sum((1-dums).*WijLPE,1)./SumWikLPE; 
    LambdaPE = SumWikHPE./SumWikLPE;
    HnPE = (Gstar.*(1-1./LambdaPE)).^(1+zstar);
    NhPE = HnPE.*N_data; NlPE = (1-HnPE).*N_data;
    MijHPE = NhPE.*(WijHPE./SumWikHPE);
    MijLPE = NlPE.*(WijLPE./SumWikLPE);
    LhPE = sum(MijHPE,2)'; LlPE = sum(MijLPE,2)';
    IhPE = sum((1-dums).*MijHPE,2)'; Il = sum((1-dums).*MijLPE,2)';
    MhPE = sum((1-dums).*MijHPE,1); MlPE = sum((1-dums).*MijLPE,1);
    LPE  = LhPE + LlPE;
    HrPE = LhPE./(LhPE+LlPE);
    NiihPE = NhPE.*((wH.^(1/mu))./SumWikHPE);
    NiilPE = NlPE.*((wL.^(1/mu))./SumWikLPE);
    mijPE = (MijHPE+MijLPE)./(NhPE+NlPE);
    mPE = sum((1-dums).*mijPE,1); % Outmigration rate New 
    % erreor computation 
    errora  = abs(LhPE - LhPEtemp); 
    errorb  = abs(LlPE - LlPEtemp); 
    error=sum(errora(:)+errorb(:));

    fprintf('%6g %6g \n',count, error);

    % assign new value / update the guess
    LhPEtemp =corr*LhPEtemp + (1-corr)*LhPE;
    LlPEtemp =corr*LlPEtemp + (1-corr)*LlPE;

    HrPE = LhPEtemp./(LhPEtemp+LlPEtemp);
    LiTPE =(thetHPE.*HrPE.^psi + thetLPE.*(1-HrPE).^psi).^(1/psi);

 end 

 %% Development accounting 
    ywdPE =TFPPE.*LiTPE.*(1-tauPE);
    NIR1_PE  = ywd./ywdPE -1; % Human capital channel
    
    SRN = HrPE./(1-HrPE);
    RPE = SRN./SRL;
    R_ThetPE = thetHPE./thetL;
    thetHPE =(R_ThetPE.*RPE.^kappa)./(1+(R_ThetPE.*RPE.^kappa));
    thetLPE =  1-thetHPE;
    LiTPE =(thetHPE.*HrPE.^psi + thetLPE.*(1-HrPE).^psi).^(1/psi);
    TFPPE = TFPPE.*RPE.^epsilon;
    ywdPE =TFPPE.*LiTPE.*(1-tauPE);
    NIR2_PE = ywd./ywdPE -1; % Technological externality 

    denomi = (thetHPE.*Hr.^psi +thetLPE.*(1-Hr).^psi).^(1/psi); % technological externality only
    numera = (thetHPE.*HrPE.^psi + thetLPE.*((1-HrPE).^psi)).^(1/psi);
    TEC_PE  = (TFP./TFPPE - 1) +(numera./denomi -1); 

    TFPPE   = TFPPE.*(((mbar+mPE)./(mbar+ mi)).^rho); % Diaspora externalities
    ywdPE   = TFPPE.*LiTPE.*(1-tauPE);
    NIR3_PE = ywd./ywdPE -1;

    TFPPE = TFPPE.*((mbar+mPE)./(mbar+ mi)).^rho; % Diaspora externality only
    DIAS_PE  = TFP./TFPPE - 1; 
 
    tauPE = GovCons_data.*LPEtemp.^(-eta) + EduExp_data./(1-mPE); % Account for fiscal externalities (tax and size)
    ywdPE= TFPPE.*LiTPE.*(1-tauPE);   
    FIS_PE  = (tauPE - tau1)./(1-tauPE); % fiscal externality only
    NIR4_PE = NIR3_PE + FIS_PE; 


    NQ =LiT ; % Account for market size effects
    NQPE  = NQ.*((1-mPE)./(1-mi));  LossNQPE =(NQ - NQPE)./NQPE;
    NIR5_PE  = NIR4_PE + 0.7*(LossNQPE./(lambda1 - 1));
    MKT_PE  = 0.7*(LossNQPE./(lambda1 - 1)); % market externality only

   
    MijPE = (1-dums).*MijHPE+(1-dums).*MijLPE;
    RijPE = Rij_pm.*MijPE;
    rPE = sum(RijPE,1)./LPE;  
    ywd_PE= ywdPE.*(1-MKT_PE);
 
    REM_PE  = (r_data - rPE)./(ywdPE+rPE); % remittances externality only
    NIR6_PE  = NIR5_PE + REM_PE;

     % Income per natural
    wijPE = dums.*(wLPE.*(1-tauPE).*NiilNew + wHPE.*(1-tauPE).*NiihNew) + (1-dums).*(wL'.*(1-tau1').*MijLPE + wH'.*(1-tau1').*MijHPE);
    ypndPE = sum(wijPE,1)./N_data; 
    dwnPE   = ypnd./ypndPE-1; % relative deviation in income per natural due to migration 
    NIR7_PE= NIR5_PE+dwnPE; % remittances are factored out as already part of migrants' income abroad
    NIR7_PE = min(.9,NIR7_PE);


%--------------------------------------------------------------------------
%   CHANNELS
%--------------------------------------------------------------------------

c_HumCap = NIR1_PE;
c_TecExt = NIR2_PE - NIR1_PE;
c_DiaExt = NIR3_PE - NIR2_PE;
c_FisExt = NIR4_PE - NIR3_PE;
c_MktExt = NIR5_PE - NIR4_PE;
c_RemEff = NIR6_PE - NIR5_PE;
c_NatEff = NIR7_PE - NIR6_PE;

cc_HumCap = NIR1_PE;
cc_TecExt = TEC_PE;
cc_DiaExt = DIAS_PE;
cc_FisExt = FIS_PE;
cc_MktExt = MKT_PE;
cc_RemEff = REM_PE;
cc_NatEff = NIR7_PE - NIR6_PE;

c_channels =[c_HumCap',c_TecExt',c_DiaExt',c_FisExt',c_MktExt',c_RemEff',c_NatEff',...
             cc_HumCap',cc_TecExt',cc_DiaExt',cc_FisExt',cc_MktExt',cc_RemEff',cc_NatEff',...
             NIR6_PE',NIR7_PE',ywd_nomig',ywd1'];
cols = ["HumCapPE","TecExtPE","DiaExtPE","FisExtPE","MktExtPE","RemEffPE","NatEffPE",...
        "HumCapPE2","TecExtPE2","DiaExtPE2","FisExtPE2","MktExtPE2","RemEffPE2","NatEffPE2",...
        "NIR6_PEb","NIR7_PEb","yPE_b","ywdPE"];
cc_channels_PE =[cols;c_channels];
