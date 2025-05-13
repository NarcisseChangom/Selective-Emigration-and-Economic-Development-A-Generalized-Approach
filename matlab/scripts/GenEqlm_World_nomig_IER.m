

%-------------------------------------------------------------------------%
%          WORLD GENERAL EQUILIBRIUM: (ENDO. WAGES  ORIG. & DEST)         %
%-------------------------------------------------------------------------%

% No selection 
global sigma psi mu kappa epsilon eta lambda1 rho mbar zstar thetHGE thetLGE
global N_data L_data Gstar HrGE HnGE LambdaGE LiTGE dums TFPGE tauGE 

%------------------------
% A: Equilibrium replication  
%------------------------ 
cijHGE = (1-dums).*ones(o,d); % Shock  
cijLGE = (1-dums).*ones(o,d); % Shock  
LiTGE = LiT; thetHGE = thetH; thetLGE = 1-thetHGE; TFPGE = TFP; tauGE = tau1;
LhGEtemp = Lh; LlGEtemp = Ll;


corr=0.9; maxerror =1.0e-5; maxit = 5000; error=2; count = 0; 
 while (error > maxerror) && (count < maxit)
      count = count + 1;
       error=0; errora=0;  errorb=0;

    % Skill specific endogeneous wages 
    LGEtemp = LhGEtemp+LlGEtemp;
    HrGE = LhGEtemp./LGEtemp;
    wHGE = thetHGE.*TFPGE.*(LiTGE./HrGE).^(1/sigma);
    wLGE = thetLGE.*TFPGE.*(LiTGE./(1-HrGE)).^(1/sigma);

    % Corridor specific utility and optimal re-allocation choice 
    WijHGE = (wHGE'.*(1-cijHGE)).^(1/mu);  
    WijLGE = (wLGE'.*(1-cijLGE)).^(1/mu); 

    SumWikHGE = wHGE.^(1/mu) + sum((1-dums).*WijHGE,1); % endogeneous wages at origin and exogenous at destination 
    SumWikLGE = wLGE.^(1/mu) + sum((1-dums).*WijLGE,1);

    mHGE = sum((1-dums).*WijHGE,1)./SumWikHGE;
    mLGE = sum((1-dums).*WijLGE,1)./SumWikLGE; 
    LambdaGE = SumWikHGE./SumWikLGE;
    HnGE = (Gstar.*(1-1./LambdaGE)).^(1+zstar);
    NhGE = HnGE.*N_data; NlGE = (1-HnGE).*N_data;
    MijHGE = NhGE.*(WijHGE./SumWikHGE);
    MijLGE = NlGE.*(WijLGE./SumWikLGE);
    LhGE = sum(MijHGE,2)'; LlGE = sum(MijLGE,2)';
    IhGE = sum((1-dums).*MijHGE,2)'; Il = sum((1-dums).*MijLGE,2)';
    MhGE = sum((1-dums).*MijHGE,1); MlPE = sum((1-dums).*MijLGE,1);
    LGE  = LhGE + LlGE;
    HrGE = LhGE./(LhGE+LlGE);
    NiihGE = NhGE.*((wH.^(1/mu))./SumWikHGE);
    NiilGE = NlGE.*((wL.^(1/mu))./SumWikLGE);

    mijGE = (MijHGE+MijLGE)./(NhGE+NlGE);
    mGE = sum((1-dums).*mijGE,1); % Outmigration rate New 

    % erreor computation 
    errora  = abs(LhGE - LhGEtemp); 
    errorb  = abs(LlGE - LlGEtemp); 
    error=sum(errora(:)+errorb(:));

    fprintf('%6g %6g \n',count, error);

    % assign new value / update the guess
    LhGEtemp =corr*LhGEtemp + (1-corr)*LhGE;
    LlGEtemp =corr*LlGEtemp + (1-corr)*LlGE;

    HrGE = LhGEtemp./(LhGEtemp+LlGEtemp);
    LiTGE =(thetHGE.*HrGE.^psi + thetLGE.*(1-HrGE).^psi).^(1/psi);

 end 

 %% Development accounting 
    ywdGE =TFPGE.*LiTGE.*(1-tauGE);
    NIR1_GE  = ywd./ywdGE -1; % Human capital channel
    
    SRNGE = HrGE./(1-HrGE);
    RGE = SRNGE./SRL;
    R_ThetGE = thetHGE./thetL;
    thetHGE =(R_ThetGE.*RGE.^kappa)./(1+(R_ThetGE.*RGE.^kappa));
    thetLGE =  1-thetHGE;
    LiTGE =(thetHGE.*HrGE.^psi + thetLGE.*(1-HrGE).^psi).^(1/psi);
    TFPGE = TFPGE.*RGE.^epsilon;
    ywdGE =TFPGE.*LiTGE.*(1-tauGE);
    NIR2_GE = ywd./ywdGE -1; % Technological externality 

    denomi = (thetHGE.*Hr.^psi +thetLGE.*(1-Hr).^psi).^(1/psi); % technological externality only
    numera = (thetHGE.*HrGE.^psi + thetLGE.*((1-HrGE).^psi)).^(1/psi);
    TEC_GE  = (TFP./TFPGE - 1) +(numera./denomi -1); 

    TFPGE   = TFPGE.*(((mbar+mGE)./(mbar+ mi)).^rho); % Diaspora externalities
    ywdGE   = TFPGE.*LiTGE.*(1-tauGE);
    NIR3_GE = ywd./ywdGE -1;

    TFPGE = TFPGE.*((mbar+mGE)./(mbar+ mi)).^rho; % Diaspora externality only
    DIAS_GE  = TFP./TFPGE - 1; 
 
    tauGE = GovCons_data.*LGEtemp.^(-eta) + EduExp_data./(1-mGE); % Account for fiscal externalities (tax and size)
    ywdGE= TFPGE.*LiTGE.*(1-tauGE);   
    FIS_GE  = (tauGE - tau1)./(1-tauGE); % fiscal externality only
    NIR4_GE = NIR3_GE + FIS_GE; 

    NQ =LiT ; % Account for market size effects
    NQGE  = NQ.*((1-mGE)./(1-mi));  LossNQGE =(NQ - NQGE)./NQGE;
    NIR5_GE  = NIR4_GE + 0.7*(LossNQGE./(lambda1 - 1));
    MKT_GE  = 0.7*(LossNQGE./(lambda1 - 1)); % market externality only
 
    %% Remittances   
    MijGE = (1-dums).*MijHGE+(1-dums).*MijLGE;
    RijGE = Rij_pm.*MijGE; % remittances per migrant at corridor level is assumed exogeneous (we do not model remitting decision of migrants)
    rGE = sum(RijGE,1)./LGE;  
    ywd_GE= ywdGE.*(1-MKT_GE);

    NIR6_GE  = NIR5_GE + (r_data - rGE)./(ywdGE+rGE);  %   Account for remittances
    REM_GE  = (r_data - rGE)./(ywdGE+rGE); % remittances externality only

     % Income per natural
    wijGE = dums.*(wLGE.*(1-tauGE).*NiilGE + wHGE.*(1-tauGE).*NiihGE) + (1-dums).*(wLGE'.*(1-tauGE').*MijLGE + wHGE'.*(1-tauGE').*MijHGE);
    ypndGE = sum(wijGE,1)./N_data; 
    dwnGE   = ypnd./ypndGE-1; % relative deviation in income per natural due to migration 
    NIR7_GE = NIR5_GE+dwnGE; % remittances are factored out as already part of migrants' income abroad
    NIR7_GE = min(.9,NIR7_GE);


  %% Income inequality (Theil index) 
  wHGEd = wHGE.*(1-tauGE);
  wLGEd = wLGE.*(1-tauGE);
  wGEd  = HrGE.*wHGEd + (1-HrGE).*wLGEd;

  
  % ConstantPop
 wdbarGE = sum(wGEd.*LGE)./sum(LGE); % World average disposable income 
 TheilGE = (sum(wHGEd.*LhGE.*log(wHGEd./wdbarGE)) + sum(wLGEd.*LlGE.*log(wLGEd./wdbarGE)))./(wdbarGE.*sum(LGE)); % Theil index 
 TheilAGE = (sum(wGEd.*LGE.*log(wGEd./wdbarGE)))./(wdbarGE.*sum(LGE)); % Theil index across countries 
 InnerTGE = (((wHGEd.*LhGE).*log(wHGEd./wGEd)) + ((wLGEd.*LlGE).*log(wLGEd./wGEd)))./(wGEd.*LGE); % inner component of the across 
 TheilWGE = sum((wGEd.*LGE.*InnerTGE))./(wdbarGE.*sum(LGE)); % Theil index across countries 

 % NewPop
 wdbarGEc = sum(wGEd.*L_data)./sum(L_data); 
 TheilGEc = (sum(wHGEd.*Lh_data.*log(wHGEd./wdbarGEc)) + sum(wLGEd.*Ll_data.*log(wLGEd./wdbarGEc)))./(wdbarGEc.*sum(L_data)); 
 TheilAGEc = (sum(wGEd.*L_data.*log(wGEd./wdbarGEc)))./(wdbarGEc.*sum(L_data)); 
 InnerTGEc = (((wHGEd.*Lh_data).*log(wHGEd./wGEd)) + ((wLGEd.*Ll_data).*log(wLGEd./wGEd)))./(wGEd.*L_data); 
 TheilWGEc = sum((wGEd.*L_data.*InnerTGEc))./(wdbarGEc.*sum(L_data));  


%% Wages with spillovers (NO MIGRATION)

    QSC  = (thetH.*HrGE.^psi + thetL.*(1-HrGE).^psi).^(1/psi);  
    wHSC = thetH.*TFP.*(QSC./HrGE).^(1/sigma);  
    wLSC = thetL.*TFP.*(QSC./(1-HrGE)).^(1/sigma); 
    SRHSC= HrGE./(1-HrGE); 
    SRL_data= Hr_data./(1-Hr_data); 
    TFP_CF = TFP.*((SRHSC./SRL_data).^epsilon).*((mbar+mGE)./(mbar+mi)).^rho;  
    wHSCa  = thetH.*TFP_CF.*(QSC./HrGE).^(1/sigma);  
    wLSCa  = thetL.*TFP_CF.*(QSC./(1-HrGE)).^(1/sigma);  
    
    Ratio_b  = thetH./(1-thetH); 
    thetH_CF = (Ratio_b.*(SRHSC./SRL_data).^kappa)./(1+Ratio_b.*(SRHSC./SRL_data).^kappa);  
    thetL_CF = 1 -thetH_CF;
    QSC_CFb  = (thetH_CF.*HrGE.^psi + thetL_CF.*(1-HrGE).^psi).^(1/psi);  
    wHSCb = thetH_CF.*TFP_CF.*(QSC_CFb./HrGE).^(1/sigma); 
    wLSCb = thetL_CF.*TFP_CF.*(QSC_CFb./(1-HrGE)).^(1/sigma);  
    NQ    = QSC_CFb;  
    NQcf  = QSC_CFb./(1-mi);
    rem_data = r_data./ywdGE;
    remSC = rGE./ywdGE;
    LossNQ= (NQ.*(1+rem_data))./NQcf - 1;
    MS    = 0.7*LossNQ./(lambda1-1);
    wHSCc = wHSCb.*(1-MS);
    wLSCc = wLSCb.*(1-MS);     
    wHSCd = wHSCc.*(1+remSC).*(1-tauGE); % final wage HS
    wLSCd = wLSCc.*(1+remSC).*(1-tauGE); % final wage LS
    wHdd = wHd+ri;
    wLdd = wLd+ri;
    clear QSC wHSC wLSC SRHSC SRL_d SRL_data TFP_CF wHSCa...
          wLSCa Ratio_b thetH_CF thetL_CF QSC_CFb wHSCb wLSCb...
          Qcf rem_data remSC wHSCc wLSCc MS
   


%--------------------------------------------------------------------------
%   CHANNELS
%--------------------------------------------------------------------------

c_HumCap = NIR1_GE;
c_TecExt = NIR2_GE - NIR1_GE;
c_DiaExt = NIR3_GE - NIR2_GE;
c_FisExt = NIR4_GE - NIR3_GE;
c_MktExt = NIR5_GE - NIR4_GE;
c_RemEff = NIR6_GE - NIR5_GE;
c_NatEff = NIR7_GE - NIR6_GE;

cc_HumCap = NIR1_GE;
cc_TecExt = TEC_GE;
cc_DiaExt = DIAS_GE;
cc_FisExt = FIS_GE;
cc_MktExt = MKT_GE;
cc_RemEff = REM_GE;
cc_NatEff = NIR7_GE - NIR6_GE;

c_channels =[c_HumCap',c_TecExt',c_DiaExt',c_FisExt',c_MktExt',c_RemEff',c_NatEff',...
             cc_HumCap',cc_TecExt',cc_DiaExt',cc_FisExt',cc_MktExt',cc_RemEff',cc_NatEff',...
             NIR6_GE',NIR7_GE',ywd_nomig',ywd1'];
cols = ["HumCapGE","TecExtGE","DiaExtGE","FisExtGE","MktExtGE","RemEffGE","NatEffGE",...
        "HumCapGE2","TecExtGE2","DiaExtGE2","FisExtGE2","MktExtGE2","RemEffGE2","NatEffGE2",...
        "NIR6_GEb","NIR7_GEb","yGE_b","ywdGE"];
cc_channels_GE =[cols;c_channels];


