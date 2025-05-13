
%-------------------------------------------------------------------------*
%      SELECTIVE MIGRATION, HUMAN CAPITAL AND DEVELOPMENT:                *
%                   A GENERALIZED APPROACH                                *
%                   Authors:
%                       - Narcisse Cha'ngom (LISER)     
%                       - Christoph Deuster (European Union)
%                       - Frédéric Docquier (LISER) 
%                       - Joel Machado (LISER) 
%     Last update: May 07, 2025                                           *
%-------------------------------------------------------------------------*

%-------------------------------------------------------------------------%
%                  Sensitivity to  Parameter values                       %
%-------------------------------------------------------------------------%
tic 
global mu sigma psi kappa epsilon eta lambda1 rho mbar zstar_lic zstar_lmic zstar zstar_umic zstar_hic

for i=1:o
    if GNIc_data(i) <1000
        LIC(i)=1;
        else LIC(i)=0;
    end
    
    if GNIc_data(i)> 1000 && GNIc_data(i)<4000
        LMIC(i)=1;
        else LMIC(i)=0;
    end
    
     if GNIc_data(i)> 4000 && GNIc_data(i)<12500
        UMIC(i)=1;
        else UMIC(i)=0;
     end
    
     if GNIc_data(i) >12500
        HIC(i)=1;
        else HIC(i)=0;
    end
end 

%------------------------------
% Sensitivity to sigma (1.5)  
%------------------------------

% Benchmark values 
mu = 1/0.7; kappa = 0.1; epsilon =0.1; eta = 0.056;  lambda1 = 8; rho = 0.032; mbar= 0.1;
zstar_lic =5.3; zstar_lmic=3.8; zstar_umic=0.0; zstar_hic=0.0;
zstar = zstar_lic*LIC + zstar_lmic*LMIC + zstar_umic*UMIC + zstar_hic*HIC;

% Updated sigma (1.5)
sigma =1.5; psi = (sigma - 1)/sigma;
    run('matlab\scripts\Migration_costs_and_model_inversion_IER.m');
    run('matlab\scripts\InitEqlm_IER.m');
    run('matlab\scripts\Partial_Eqlm_nomig_IER.m');
    run('matlab\scripts\Partial_Eqlm_noselec_IER.m');
    run('matlab\scripts\GenEqlm_SODC_nomig_IER.m');
    run('matlab\scripts\GenEqlm_SODC_noselec_IER.m');
    run('matlab\scripts\GenEqlm_World_nomig_IER.m');
    run('matlab\scripts\GenEqlm_World_noselec_IER.m');

NetIncomeResponse = table(isoo,NIR6_NM',NIR6_NS',NIR6_PE',NIR6_PE4',NIR6_GE',NIR6_GE4',NIR7_NM',NIR7_NS',NIR7_PE',NIR7_PE4',NIR7_GE',NIR7_GE4',NIR5_NM',NIR5_PE');
NetIncomeResponse = renamevars(NetIncomeResponse,["isoo","Var2","Var3","Var4","Var5","Var6","Var7", ...
                                              "Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15"]...
                                              ,["iso","NIR6_NM","NIR6_NS","NIR6_PESC","NIR6_PESC4","NIR6_SC", ...
                                                "NIR6_SC4","NIR7_NM","NIR7_NS","NIR7_PESC","NIR7_PESC4", ...
                                                "NIR7_SC","NIR7_SC4","NIR5_NM","NIR5_PESC"]);
writetable(NetIncomeResponse, 'matlab\output\Benchmarksigma15_IER_25.csv', 'Delimiter',',' ,'QuoteStrings', true)


%------------------------------
% Sensitivity to mu  (1/0.6) 
%------------------------------

% Benchmark values 
sigma =2; psi = (sigma - 1)/sigma; kappa = 0.1; epsilon =0.1; eta = 0.056;  lambda1 = 8; rho = 0.032; mbar= 0.1;
zstar_lic =5.3; zstar_lmic=3.8; zstar_umic=0.0; zstar_hic=0.0;
zstar = zstar_lic*LIC + zstar_lmic*LMIC + zstar_umic*UMIC + zstar_hic*HIC;

% Updated mu (1/0.6)
mu = 1/0.6; 
    run('matlab\scripts\Migration_costs_and_model_inversion_IER.m');
    run('matlab\scripts\InitEqlm_IER.m');
    run('matlab\scripts\Partial_Eqlm_nomig_IER.m');
    run('matlab\scripts\Partial_Eqlm_noselec_IER.m');
    run('matlab\scripts\GenEqlm_SODC_nomig_IER.m');
    run('matlab\scripts\GenEqlm_SODC_noselec_IER.m');
    run('matlab\scripts\GenEqlm_World_nomig_IER.m');
    run('matlab\scripts\GenEqlm_World_noselec_IER.m');

NetIncomeResponse = table(isoo,NIR6_NM',NIR6_NS',NIR6_PE',NIR6_PE4',NIR6_GE',NIR6_GE4',NIR7_NM',NIR7_NS',NIR7_PE',NIR7_PE4',NIR7_GE',NIR7_GE4',NIR5_NM',NIR5_PE');
NetIncomeResponse = renamevars(NetIncomeResponse,["isoo","Var2","Var3","Var4","Var5","Var6","Var7", ...
                                              "Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15"]...
                                              ,["iso","NIR6_NM","NIR6_NS","NIR6_PESC","NIR6_PESC4","NIR6_SC", ...
                                                "NIR6_SC4","NIR7_NM","NIR7_NS","NIR7_PESC","NIR7_PESC4", ...
                                                "NIR7_SC","NIR7_SC4","NIR5_NM","NIR5_PESC"]);
writetable(NetIncomeResponse, 'matlab\output\Benchmarkmu06_IER_25.csv', 'Delimiter',',' ,'QuoteStrings', true)       



%------------------------------
% Sensitivity to mu  (0.7) 
%------------------------------

% Benchmark values 
sigma =2; psi = (sigma - 1)/sigma; kappa = 0.1; epsilon =0.1; eta = 0.056;  lambda1 = 8; rho = 0.032; mbar= 0.1;
zstar_lic =5.3; zstar_lmic=3.8; zstar_umic=0.0; zstar_hic=0.0;
zstar = zstar_lic*LIC + zstar_lmic*LMIC + zstar_umic*UMIC + zstar_hic*HIC;

% Updated mu (0.7)
mu = 0.7; 
    run('matlab\scripts\Migration_costs_and_model_inversion_IER.m');
    run('matlab\scripts\InitEqlm_IER.m');
    run('matlab\scripts\Partial_Eqlm_nomig_IER.m');
    run('matlab\scripts\Partial_Eqlm_noselec_IER.m');
    run('matlab\scripts\GenEqlm_SODC_nomig_IER.m');
    run('matlab\scripts\GenEqlm_SODC_noselec_IER.m');
    run('matlab\scripts\GenEqlm_World_nomig_IER.m');
    run('matlab\scripts\GenEqlm_World_noselec_IER.m');

NetIncomeResponse = table(isoo,NIR6_NM',NIR6_NS',NIR6_PE',NIR6_PE4',NIR6_GE',NIR6_GE4',NIR7_NM',NIR7_NS',NIR7_PE',NIR7_PE4',NIR7_GE',NIR7_GE4',NIR5_NM',NIR5_PE');
NetIncomeResponse = renamevars(NetIncomeResponse,["isoo","Var2","Var3","Var4","Var5","Var6","Var7", ...
                                              "Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15"]...
                                              ,["iso","NIR6_NM","NIR6_NS","NIR6_PESC","NIR6_PESC4","NIR6_SC", ...
                                                "NIR6_SC4","NIR7_NM","NIR7_NS","NIR7_PESC","NIR7_PESC4", ...
                                                "NIR7_SC","NIR7_SC4","NIR5_NM","NIR5_PESC"]);
writetable(NetIncomeResponse, 'matlab\output\Benchmarkmu142_IER_25.csv', 'Delimiter',',' ,'QuoteStrings', true)      



%------------------------------
% Sensitivity to z star  (SR) 
%------------------------------

% Benchmark values 
mu = 1/0.7; sigma =2; psi = (sigma - 1)/sigma; kappa = 0.1; epsilon =0.1; eta = 0.056;  lambda1 = 8; rho = 0.032; mbar= 0.1;

% Updated zstar (SR)
 zstar_lic =1.7; zstar_lmic=0.9; zstar_umic=0.0; zstar_hic=0.0;
 zstar = zstar_lic*LIC + zstar_lmic*LMIC + zstar_umic*UMIC + zstar_hic*HIC;
    run('matlab\scripts\Migration_costs_and_model_inversion_IER.m');
    run('matlab\scripts\InitEqlm_IER.m');
    run('matlab\scripts\Partial_Eqlm_nomig_IER.m');
    run('matlab\scripts\Partial_Eqlm_noselec_IER.m');
    run('matlab\scripts\GenEqlm_SODC_nomig_IER.m');
    run('matlab\scripts\GenEqlm_SODC_noselec_IER.m');
    run('matlab\scripts\GenEqlm_World_nomig_IER.m');
    run('matlab\scripts\GenEqlm_World_noselec_IER.m');

NetIncomeResponse = table(isoo,NIR6_NM',NIR6_NS',NIR6_PE',NIR6_PE4',NIR6_GE',NIR6_GE4',NIR7_NM',NIR7_NS',NIR7_PE',NIR7_PE4',NIR7_GE',NIR7_GE4',NIR5_NM',NIR5_PE');
NetIncomeResponse = renamevars(NetIncomeResponse,["isoo","Var2","Var3","Var4","Var5","Var6","Var7", ...
                                              "Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15"]...
                                              ,["iso","NIR6_NM","NIR6_NS","NIR6_PESC","NIR6_PESC4","NIR6_SC", ...
                                                "NIR6_SC4","NIR7_NM","NIR7_NS","NIR7_PESC","NIR7_PESC4", ...
                                                "NIR7_SC","NIR7_SC4","NIR5_NM","NIR5_PESC"]);
writetable(NetIncomeResponse, 'matlab\output\BenchmarkzSR_IER_25.csv', 'Delimiter',',' ,'QuoteStrings', true)    

%------------------------------
% HumanCap results export
%------------------------------

HumanCap = table(isoo,Lambda',LambdaNM',LambdaPE',LambdaGE',Hr_data',Hn_data',HrNM',HrPE',HrGE',HnNM',HnPE',HnGE');
HumanCap = renamevars(HumanCap,["isoo","Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13"]...
                                              ,["iso","Lambdap","LambdaNMp","LambdaPESCp","LambdaSCp","Hr","Hn","HrNMp","HrPESCp","HrSCp","HnNM","HnPE","HnGE"]);
writetable(HumanCap, 'matlab\output\HumanCap_SR_25.csv', 'Delimiter',',' ,'QuoteStrings', true)


%------------------------------
% Pessimistic/conservative parameter values
%------------------------------

%Benchmark values 
mu = 1/0.7; sigma =2; psi = (sigma - 1)/sigma; mbar= 0.1;

% Conservative parameter values 
kappa = 0.2; epsilon =0.2; eta = 0.112;  lambda1 = 4; rho = 0.016;
zstar_lic =1.7; zstar_lmic=0.9; zstar_umic=0.0; zstar_hic=0.0;
zstar = zstar_lic*LIC + zstar_lmic*LMIC + zstar_umic*UMIC + zstar_hic*HIC;


    run('matlab\scripts\Migration_costs_and_model_inversion_IER.m');
    run('matlab\scripts\InitEqlm_IER.m');
    run('matlab\scripts\Partial_Eqlm_nomig_IER.m');

    % Averages 
        AVG_NIR6_NMp = [NIR6_NMworld;NIR6_NMLIC;NIR6_NMLMIC;NIR6_NMUMIC;NIR6_NMHIC];
        AVG_NIR7_NMp = [NIR7_NMworld;NIR7_NMLIC;NIR7_NMLMIC;NIR7_NMUMIC;NIR7_NMHIC];
        
    run('matlab\scripts\Partial_Eqlm_noselec_IER.m');
        AVG_NIR6_NSp = [NIR6_NSworld;NIR6_NSLIC;NIR6_NSLMIC;NIR6_NSUMIC;NIR6_NSHIC];
        AVG_NIR7_NSp = [NIR7_NSworld;NIR7_NSLIC;NIR7_NSLMIC;NIR7_NSUMIC;NIR7_NSHIC];
        AVG_NIRp = table(AVG_NIR6_NMp,AVG_NIR7_NMp,AVG_NIR6_NSp,AVG_NIR7_NSp);


    run('matlab\scripts\GenEqlm_SODC_nomig_IER.m');
    run('matlab\scripts\GenEqlm_SODC_noselec_IER.m');
    run('matlab\scripts\GenEqlm_World_nomig_IER.m');
    run('matlab\scripts\GenEqlm_World_noselec_IER.m');
    
NetIncomeResponse = table(isoo,NIR6_NM',NIR6_NS',NIR6_PE',NIR6_PE4',NIR6_GE',NIR6_GE4',NIR7_NM',NIR7_NS',NIR7_PE',NIR7_PE4',NIR7_GE',NIR7_GE4',NIR5_NM',NIR5_PE');
NetIncomeResponse = renamevars(NetIncomeResponse,["isoo","Var2","Var3","Var4","Var5","Var6","Var7", ...
                                              "Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15"]...
                                              ,["iso","NIR6_NM","NIR6_NS","NIR6_PESC","NIR6_PESC4","NIR6_SC", ...
                                                "NIR6_SC4","NIR7_NM","NIR7_NS","NIR7_PESC","NIR7_PESC4", ...
                                                "NIR7_SC","NIR7_SC4","NIR5_NM","NIR5_PESC"]);
writetable(NetIncomeResponse, 'matlab\output\Pessimistic_IER_25.csv', 'Delimiter',',' ,'QuoteStrings', true)    


%% World distribution of income (conservative)    

WAGESp = table(isoo,wHdd',wLdd',wHSCd',wLSCd',LhGE',LlGE',Lh_data',Ll_data');

WAGESp = renamevars(WAGESp,["isoo","Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9"]...
                          ,["iso","wHdd","wLdd","wHSCd","wLSCd","LhSC","LlSC","Lh","Ll"]);

% Export to CSV 
writetable(WAGESp, 'matlab\output\WAGESp_25.csv', 'Delimiter',',' ,'QuoteStrings', true)

%% Channels 
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
cols = ["HumCap_p","TecExt_p","DiaExt_p","FisExt_p","MktExt_p","RemEff_p","NatEff_p","HumCap2_p","TecExt2_p","DiaExt2_p","FisExt2_p","MktExt2_p","RemEff2_p","NatEff2_p","NIR6_NMp","NIR7_NMp","ynomig_p","ywdp"];
cc_channels =[cols;c_channels];
writematrix(cc_channels, 'matlab\output\cc_channelsNM_p_25.csv', 'Delimiter',',' ,'QuoteStrings', true)

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
cols = ["HumCap_p","TecExt_p","DiaExt_p","FisExt_p","MktExt_p","RemEff_p","NatEff_p","HumCap2_p","TecExt2_p","DiaExt2_p","FisExt2_p","MktExt2_p","RemEff2_p","NatEff2_p","NIR6_NSp","NIR7_NSp","ynoselec_p","ywdp"];
cc_channels_NS =[cols;c_channels];
writematrix(cc_channels_NS, 'matlab\output\cc_channelsNS_p_25.csv', 'Delimiter',',' ,'QuoteStrings', true)


%------------------------
% Inequality results    
%------------------------ 
THEIL_Obs =[Theil;TheilA;TheilW];
THEIL_NMc =[TheilGE;TheilAGE;TheilWGE];
THEIL_NM =[TheilGEc;TheilAGEc;TheilWGEc];
THEIL_NSc =[TheilGE4;TheilAGE4;TheilWGE4];
THEIL_NS =[TheilGE4c;TheilAGE4c;TheilWGE4c];

THEIL_Pess = table(THEIL_Obs,THEIL_NMc,THEIL_NM,THEIL_NSc,THEIL_NS);
writetable(THEIL_Pess, 'matlab\output\THEIL_Pess_25.csv', 'Delimiter',',' ,'QuoteStrings', true)

THEIL_Pess = renamevars(THEIL_Pess,["THEIL_Obs","THEIL_NMc","THEIL_NM","THEIL_NSc","THEIL_NS"]...
                                  ,["THEIL_Obs_p","THEIL_NMc_p","THEIL_NM_p","THEIL_NSc_p","THEIL_NS_p"]);

THEIL_IER = [THEIL_Bench,THEIL_Pess];
THEIL_IER.THEIL_Obs_p = [];
writetable(THEIL_IER, 'matlab\output\THEIL_IER_25.csv', 'Delimiter',',' ,'QuoteStrings', true)

%------------------------
% Averages (NIR) results (Bench+Conservative)   
%------------------------ 
AVG_NIR = [AVG_NIRb,AVG_NIRp];

writetable(AVG_NIR, 'matlab\output\AVG_NIR_25.csv', 'Delimiter',',' ,'QuoteStrings', true)

toc 