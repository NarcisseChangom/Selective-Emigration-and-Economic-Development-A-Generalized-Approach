

%-------------------------------------------------------------------------%
%          WORLD GENERAL EQUILIBRIUM: (ENDO. WAGES  ORIG. & DEST)         %
%-------------------------------------------------------------------------%

% No selection 
global sigma psi mu kappa epsilon eta lambda1 rho mbar zstar thetHGE4 thetLGE4
global N_data L_data Gstar HrGE4 HnGE4 LambdaGE4 LiTGE4 dums TFPGE4 tauGE4 

%------------------------
% A: Equilibrium replication  
%------------------------ 
cijHGE4 = (1-dums).*cijH_NS; cijLGE4 =  (1-dums).*cijL_NS; 
LiTGE4 = LiT; thetHGE4 = thetH; thetLGE4 = 1-thetHGE4; TFPGE4 = TFP; tauGE4 = tau1;
LhGE4temp = Lh; LlGE4temp = Ll;


corr=0.9; maxerror =1.0e-5; maxit = 5000; error=2; count = 0; 
 while (error > maxerror) && (count < maxit)
      count = count + 1;
       error=0; errora=0;  errorb=0;

    % Skill specific endogeneous wages 
    LGE4temp = LhGE4temp+LlGE4temp;
    HrGE4 = LhGE4temp./LGE4temp;
    wHGE4 = thetHGE4.*TFPGE4.*(LiTGE4./HrGE4).^(1/sigma);
    wLGE4 = thetLGE4.*TFPGE4.*(LiTGE4./(1-HrGE4)).^(1/sigma);

    % Corridor specific utility and optimal re-allocation choice 
    WijHGE4 = (wHGE4'.*(1-cijHGE4)).^(1/mu);  
    WijLGE4 = (wLGE4'.*(1-cijLGE4)).^(1/mu); 

    SumWikHGE4 = wHGE4.^(1/mu) + sum((1-dums).*WijHGE4,1); % endogeneous wages at origin and exogenous at destination 
    SumWikLGE4 = wLGE4.^(1/mu) + sum((1-dums).*WijLGE4,1);

    mHGE4 = sum((1-dums).*WijHGE4,1)./SumWikHGE4;
    mLGE4 = sum((1-dums).*WijLGE4,1)./SumWikLGE4; 
    LambdaGE4 = SumWikHGE4./SumWikLGE4;
    HnGE4 = (Gstar.*(1-1./LambdaGE4)).^(1+zstar);
    NhGE4 = HnGE4.*N_data; NlGE4 = (1-HnGE4).*N_data;
    MijHGE4 = NhGE4.*(WijHGE4./SumWikHGE4);
    MijLGE4 = NlGE4.*(WijLGE4./SumWikLGE4);
    LhGE4 = sum(MijHGE4,2)'; LlGE4 = sum(MijLGE4,2)';
    IhGE4 = sum((1-dums).*MijHGE4,2)'; Il = sum((1-dums).*MijLGE4,2)';
    MhGE4 = sum((1-dums).*MijHGE4,1); MlPE = sum((1-dums).*MijLGE4,1);
    LGE4  = LhGE4 + LlGE4;
    HrGE4 = LhGE4./(LhGE4+LlGE4);
    NiihGE4 = NhGE4.*((wH.^(1/mu))./SumWikHGE4);
    NiilGE4 = NlGE4.*((wL.^(1/mu))./SumWikLGE4);

    mijGE4 = (MijHGE4+MijLGE4)./(NhGE4+NlGE4);
    mGE4 = sum((1-dums).*mijGE4,1); % Outmigration rate New 

    % erreor computation 
    errora  = abs(LhGE4 - LhGE4temp); 
    errorb  = abs(LlGE4 - LlGE4temp); 
    error=sum(errora(:)+errorb(:));

    fprintf('%6g %6g \n',count, error);

    % assign new value / update the guess
    LhGE4temp =corr*LhGE4temp + (1-corr)*LhGE4;
    LlGE4temp =corr*LlGE4temp + (1-corr)*LlGE4;

    HrGE4 = LhGE4temp./(LhGE4temp+LlGE4temp);
    LiTGE4 =(thetHGE4.*HrGE4.^psi + thetLGE4.*(1-HrGE4).^psi).^(1/psi);

 end 

 %% Development accounting 
    ywdGE4 =TFPGE4.*LiTGE4.*(1-tauGE4);
    NIR1_GE4  = ywd./ywdGE4 -1; % Human capital channel
    
    SRNGE4 = HrGE4./(1-HrGE4);
    RGE4 = SRNGE4./SRL;
    R_ThetGE4 = thetHGE4./thetLGE4;
    thetHGE4 =(R_ThetGE4.*RGE4.^kappa)./(1+(R_ThetGE4.*RGE4.^kappa));
    thetLGE4 =  1-thetHGE4;
    LiTGE4 =(thetHGE4.*HrGE4.^psi + thetLGE4.*(1-HrGE4).^psi).^(1/psi);
    TFPGE4 = TFPGE4.*RGE4.^epsilon;
    ywdGE4 =TFPGE4.*LiTGE4.*(1-tauGE4);
    NIR2_GE4 = ywd./ywdGE4 -1; % Technological externality 

    denomi = (thetHGE4.*Hr.^psi +thetLGE4.*(1-Hr).^psi).^(1/psi); % technological externality only
    numera = (thetHGE4.*HrGE4.^psi + thetLGE4.*((1-HrGE4).^psi)).^(1/psi);
    TEC_GE4  = (TFP./TFPGE4 - 1) +(numera./denomi -1); 

    TFPGE4   = TFPGE4.*(((mbar+mGE4)./(mbar+ mi)).^rho); % Diaspora externalities
    ywdGE4   = TFPGE4.*LiTGE4.*(1-tauGE4);
    NIR3_GE4 = ywd./ywdGE4 -1;

    TFPGE4 = TFPGE4.*((mbar+mGE4)./(mbar+ mi)).^rho; % Diaspora externality only
    DIAS_GE  = TFP./TFPGE4 - 1; 
 
    tauGE4 = GovCons_data.*LGE4temp.^(-eta) + EduExp_data./(1-mGE4); % Account for fiscal externalities (tax and size)
    ywdGE4= TFPGE4.*LiTGE4.*(1-tauGE4);   
    FIS_GE4  = (tauGE4 - tau1)./(1-tauGE4); % fiscal externality only
    NIR4_GE4    = ywd./ywdGE4 - 1;

    NQ =LiT ; % Account for market size effects
    NQGE4  = NQ.*((1-mGE4)./(1-mi));
    LossNQGE =(NQ - NQGE4)./NQGE4;
    NIR5_GE4  = NIR4_GE4 + 0.7*(LossNQGE./(lambda1 - 1));
    MKT_GE4  = 0.7*(LossNQGE./(lambda1 - 1)); % market externality only
 
    %% Remittances   
    MijGE4 = (1-dums).*MijHGE4+(1-dums).*MijLGE4;
    RijGE4 = Rij_pm.*MijGE4; % remittances per migrant at corridor level is assumed exogeneous (we do not model remitting decision of migrants)
    rGE4 = sum(RijGE4,1)./LGE4;  
    ywd_GE4= ywdGE4.*(1-MKT_GE4);
    ywdGE4 = ywdGE4+rGE4;
    REM_GE4  = (r_data - rGE4)./(ywdGE4+rGE4); % remittances externality only
    NIR6_GE4  = ywd2./ywdGE4-1 + MKT_GE4;

     % Income per natural
    wijGE4 = dums.*(wLGE4.*(1-tauGE4).*NlGE4 + wHGE4.*(1-tauGE4).*NhGE4);% + (1-dums).*(wLGE4'.*(1-tauGE4').*MijLGE4 + wHGE4'.*(1-tauGE4').*MijHGE4);
    ypndGE4 = sum(wijGE4,1)./N_data; 
    dwnGE4   = ypnd./ypndGE4-1; % relative deviation in income per natural due to migration 
    NIR7_GE4 = NIR5_GE+dwnGE4; % remittances are factored out as already part of migrants' income abroad
    NIR7_GE4 = min(.9,NIR7_GE4);


 %% Income inequality (Theil index) 
 wHGE4d = wHGE4.*(1-tauGE4);
 wLGE4d = wLGE4.*(1-tauGE4);
 wGE4d  = HrGE4.*wHGE4d + (1-HrGE4).*wLGE4d;

    % NewPop
 wdbarGE4 = sum(wGE4d.*LGE4)./sum(LGE4); % World average disposable income 
 TheilGE4 = (sum(wHGE4d.*LhGE4.*log(wHGE4d./wdbarGE4)) + sum(wLGE4d.*LlGE4.*log(wLGE4d./wdbarGE4)))./(wdbarGE4.*sum(LGE4)); % Theil index 
 TheilAGE4 = (sum(wGE4d.*LGE4.*log(wGE4d./wdbarGE4)))./(wdbarGE4.*sum(LGE4)); % Theil index across countries 
 InnerTGE4 = (((wHGE4d.*LhGE4).*log(wHGE4d./wGE4d)) + ((wLGE4d.*LlGE4).*log(wLGE4d./wGE4d)))./(wGE4d.*LGE4); % inner component of the across 
 TheilWGE4 = sum((wGE4d.*LGE4.*InnerTGE4))./(wdbarGE4.*sum(LGE4)); % Theil index across countries 


 % ConstantPop
 wdbarGE4c = sum(wGE4d.*L_data)./sum(L_data); % World average disposable income
 TheilGE4c = (sum(wHGE4d.*Lh_data.*log(wHGE4d./wdbarGE4c)) + sum(wLGE4d.*Ll_data.*log(wLGE4d./wdbarGE4c)))./(wdbarGE4c.*sum(L_data)); % Theil index 
 TheilAGE4c = (sum(wGE4d.*L_data.*log(wGE4d./wdbarGE4c)))./(wdbarGE4c.*sum(L_data)); % Theil index across countries 
 InnerTGE4c = (((wHGE4d.*Lh_data).*log(wHGE4d./wGE4d)) + ((wLGE4d.*Ll_data).*log(wLGE4d./wGE4d)))./(wGE4d.*L_data); % inner component of the across 
 TheilWGE4c = sum((wGE4d.*L_data.*InnerTGE4c))./(wdbarGE4c.*sum(L_data)); % Theil index across countries 


 %% Wages with spillovers (NO MIGRATION)

    QSC4  = (thetH.*HrGE4.^psi + thetL.*(1-HrGE4).^psi).^(1/psi);  
    wHSC4 = thetH.*TFP.*(QSC4./HrGE4).^(1/sigma);  
    wLSC4 = thetL.*TFP.*(QSC4./(1-HrGE4)).^(1/sigma); 
    SRHSC4= HrGE4./(1-HrGE4); 
    SRL_data= Hr_data./(1-Hr_data); 
    TFP_CF4 = TFP.*((SRHSC4./SRL_data).^epsilon).*((mbar+mGE4)./(mbar+mi)).^rho;  
    wHSC4a  = thetH.*TFP_CF4.*(QSC4./HrGE4).^(1/sigma);  
    wLSC4a  = thetL.*TFP_CF4.*(QSC4./(1-HrGE4)).^(1/sigma);  
    
    Ratio_b4  = thetH./(1-thetH); 
    thetH_CF4 = (Ratio_b4.*(SRHSC4./SRL_data).^kappa)./(1+Ratio_b4.*(SRHSC4./SRL_data).^kappa);  
    thetL_CF4 = 1 -thetH_CF4;
    QSC_CFb4  = (thetH_CF4.*HrGE4.^psi + thetL_CF4.*(1-HrGE4).^psi).^(1/psi);  
    wHSC4b = thetH_CF4.*TFP_CF4.*(QSC_CFb4./HrGE4).^(1/sigma); 
    wLSC4b = thetL_CF4.*TFP_CF4.*(QSC_CFb4./(1-HrGE4)).^(1/sigma);  
    NQ    = QSC_CFb4;  
    NQcf  = QSC_CFb4.*((1-mGE4)./(1-mi));
    rem_data = r_data./ywdGE4;
    remSC4 = rGE4./ywdGE4;
    LossNQ4= (QSC_CFb4 - NQcf)./NQGE4;
    MS4    = 0.7*LossNQ4./(lambda1-1);
    wHSC4c = wHSC4b.*(1-MS4);
    wLSC4c = wLSC4b.*(1-MS4);     
    wHSC4d = wHSC4c.*(1+remSC4).*(1-tauGE); % final wage HS
    wLSC4d = wLSC4c.*(1+remSC4).*(1-tauGE); % final wage LS

    clear QSC4 wHSC4 wLSC4 SRHSC4 SRL_d SRL_data TFP_CF4 wHSC4a...
          wLSC4a Ratio_b4 thetH_CF4 thetL_CF4 QSC_CFb4 wHSC4b wLSC4b...
          Qcf rem_data remSC4 wHSC4c wLSC4c MS4

%--------------------------------------------------------------------------
%   CHANNELS
%--------------------------------------------------------------------------

c_HumCap = NIR1_GE4;
c_TecExt = NIR2_GE4 - NIR1_GE4;
c_DiaExt = NIR3_GE4 - NIR2_GE4;
c_FisExt = NIR4_GE4 - NIR3_GE4;
c_MktExt = NIR5_GE4 - NIR4_GE4;
c_RemEff = NIR6_GE4 - NIR5_GE4;
c_NatEff = NIR7_GE4 - NIR6_GE4;

cc_HumCap = NIR1_GE4;
cc_TecExt = TEC_GE4;
cc_DiaExt = DIAS_GE;
cc_FisExt = FIS_GE4;
cc_MktExt = MKT_GE4;
cc_RemEff = REM_GE4;
cc_NatEff = NIR7_GE4 - NIR6_GE4;

c_channels =[c_HumCap',c_TecExt',c_DiaExt',c_FisExt',c_MktExt',c_RemEff',c_NatEff',...
             cc_HumCap',cc_TecExt',cc_DiaExt',cc_FisExt',cc_MktExt',cc_RemEff',cc_NatEff',...
             NIR6_GE4',NIR7_GE4',ywd_nomig',ywd1'];
cols = ["HumCapGE4","TecExtGE4","DiaExtGE4","FisExtGE4","MktExtGE4","RemEffGE4","NatEffGE4",...
        "HumCapGE24","TecExtGE24","DiaExtGE24","FisExtGE24","MktExtGE24","RemEffGE24","NatEffGE24",...
        "NIR6_GEb","NIR7_GEb","yGE_b","ywdGE"];
cc_channels_GE4 =[cols;c_channels];

