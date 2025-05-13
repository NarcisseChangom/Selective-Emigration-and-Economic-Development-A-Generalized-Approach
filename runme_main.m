
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
clear all 
tic
% Set directory 
setpref('MyApp', 'mydir', 'add directory here\Replication-ier');
mydir = getpref('MyApp', 'mydir');
cd(mydir)

 
%%%
% Load raw data 
%%%
migdta = readtable('matlab\data\Movers_IER.csv','Headerlines',0,'ReadVariableNames',true);
migdta = sortrows(migdta,{'isoo','isod'}); 

locationdta = readtable('matlab\data\Location_charac_IER.csv','Headerlines',0,'ReadVariableNames',true);
locationdta = sortrows(locationdta,{'isoo'}); 


%% Benchmark (Master program) 
run('matlab\scripts\dataprepa_IER.m'); % Load data 

run('matlab\scripts\Parameters_Bench_IER.m'); % Set parameter values 

run('matlab\scripts\Monte_carlo_like_calibration_z_IER.m'); % Identifies the optimal z by level of development 
saveas(gcf, 'main\figures\Figure_B1.png');   

run('matlab\scripts\Migration_costs_and_model_inversion_IER.m'); % Calibrate the migration costs and invert the model to recover unobservable location characteristics 
% Export to CSV 
writetable(MIGCOST, 'matlab\output\MIGCOST.csv', 'Delimiter',',' ,'QuoteStrings', true)

run('matlab\scripts\InitEqlm_IER.m'); % Prepare the initial equilibrium 


run('matlab\scripts\Partial_Eqlm_nomig_IER.m'); % Simulate a no-migration counterfactual under exogeneous wages 
    AVG_NIR6_NM = [NIR6_NMworld;NIR6_NMLIC;NIR6_NMLMIC;NIR6_NMUMIC;NIR6_NMHIC];
    AVG_NIR7_NM = [NIR7_NMworld;NIR7_NMLIC;NIR7_NMLMIC;NIR7_NMUMIC;NIR7_NMHIC];
run('matlab\scripts\Partial_Eqlm_noselec_IER.m'); % Simulate neutral selection counterfactual under exogeneous wages 
    AVG_NIR6_NS = [NIR6_NSworld;NIR6_NSLIC;NIR6_NSLMIC;NIR6_NSUMIC;NIR6_NSHIC];
    AVG_NIR7_NS = [NIR7_NSworld;NIR7_NSLIC;NIR7_NSLMIC;NIR7_NSUMIC;NIR7_NSHIC];
    AVG_NIRb = table(AVG_NIR6_NM,AVG_NIR7_NM,AVG_NIR6_NS,AVG_NIR7_NS);

run('matlab\scripts\GenEqlm_SODC_nomig_IER.m'); % Simulate a no-migration counterfactual under endogeneous wages at origin and exogeneous wages at destination (SODC hypothesis)
run('matlab\scripts\GenEqlm_SODC_noselec_IER.m'); % Simulate a neutral selection counterfactual under endogeneous wages at origin and exogeneous wages at destination (SODC hypothesis)


run('matlab\scripts\GenEqlm_World_nomig_IER.m');  % Simulate a no-migration counterfactual under endogeneous wages at origin and destination (World General Equilibrium)
run('matlab\scripts\GenEqlm_World_noselec_IER.m');  % Simulate a neutral selection counterfactual under endogeneous wages at origin and destination (World General Equilibrium)


%------------------------
% Averages (NIR) results    
%------------------------ 
% U:\CHANGOM\CDDM2021_Completefile\CHAP1\Replication file\CALIBOUT
writetable(AVG_NIRb, 'matlab\output\AVG_NIR_Bench_25.csv', 'Delimiter',',' ,'QuoteStrings', true)


%% Channels 
writematrix(cc_channels, 'matlab\output\cc_channelsNM_b_25.csv', 'Delimiter',',' ,'QuoteStrings', true)
writematrix(cc_channels_NS, 'matlab\output\cc_channelsNS_b_25.csv', 'Delimiter',',' ,'QuoteStrings', true)


%------------------------
% Benchmark results    
%------------------------ 
% NIR6 & NIR7

NetIncomeResponse = table(isoo,NIR6_NM',NIR6_NS',NIR6_PE',NIR6_PE4',NIR6_GE',NIR6_GE4',NIR7_NM',NIR7_NS',NIR7_PE',NIR7_PE4',NIR7_GE',NIR7_GE4',NIR5_NM',NIR5_PE');
NetIncomeResponse = renamevars(NetIncomeResponse,["isoo","Var2","Var3","Var4","Var5","Var6","Var7", ...
                                              "Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15"]...
                                              ,["iso","NIR6_NM","NIR6_NS","NIR6_PESC","NIR6_PESC4","NIR6_SC", ...
                                                "NIR6_SC4","NIR7_NM","NIR7_NS","NIR7_PESC","NIR7_PESC4", ...
                                                "NIR7_SC","NIR7_SC4","NIR5_NM","NIR5_PESC"]);
% Export to CSV 
writetable(NetIncomeResponse, 'matlab\output\Benchmark_IER_25.csv', 'Delimiter',',' ,'QuoteStrings', true)


%------------------------
% World distribution of income     
%------------------------
lic = LIC';
lmic = LMIC';
umic = UMIC';
hic  = HIC';
WAGES = table(isoo,mH',mL',mi',ywd',ywdGE',w',wH',wL',wHd',wLd',wHdd',wLdd'...
             ,wHGEd',wLGEd',wHSCd',wLSCd',wHGE4d',wLGE4d',wHSC4d',wLSC4d'...
             ,ywd_PE',ywdGE',ypnd',LhNM',LlNM',LNM',LhNS',LlNS',LNS',LhPE'...
             ,LlPE',LPE',LhGE',LlGE',LGE',Lh_data',Ll_data',L_data',wHGE',wLGE' ...
             ,wHGE4',wLGE4',r_data',Nh',N_data',thetH',thetHNew',TFP',TFPNew' ...
             ,LiT',LiTNew',lic,lmic,umic,hic);

WAGES = renamevars(WAGES,["isoo","Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13"...
                         ,"Var14","Var15","Var16","Var17","Var18","Var19","Var20","Var21"...
                         ,"Var22","Var23","Var24","Var25","Var26","Var27","Var28","Var29","Var30","Var31"...
                         ,"Var32","Var33","Var34","Var35","Var36","Var37","Var38","Var39","Var40","Var41" ...
                         ,"Var42","Var43","Var44","Var45","Var46","Var47","Var48","Var49","Var50","Var51","Var52"]...
                         ,["iso","mH","mL","mi","ywd","ypwdSC","w","wH","wL","wHd","wLd","wHdd","wLdd"...
                         ,"wHGEd","wLGEd","wHSCd","wLSCd","wHGE4d","wLGE4d","wHSC4d","wLSC4d"...
                         ,"ywdPESC","ywdSC","ypnd","LhNM","LlNM","LNM","LhNS","LlNS","LNS","LhPESC"...
                         ,"LlPESC","LPESC","LhSC","LlSC","LSC","Lh","Ll","L","wHGE","wLGE","wHGE4","wLGE4" ...
                         ,"rpc","Nh","N","thetH","thetHNew","TFP","TFPNew","LiT","LiTNew"]);


% Export to CSV 
writetable(WAGES, 'matlab\output\WAGES_25.csv', 'Delimiter',',' ,'QuoteStrings', true)


% Access to public education Gstar 
G = Gstar';
gnic = GNIc_data';
PuEduc = table(isoo,G,gnic,lic,lmic,umic,hic);
writetable(PuEduc, 'matlab\output\GstarValidation_25.csv', 'Delimiter',',' ,'QuoteStrings', true)

%------------------------
% HumanCap results 
%------------------------ 

HumanCap = table(isoo,Lambda',LambdaNM',LambdaPE',LambdaGE',Hr_data'...
                ,Hn_data',HrNM',HrPE',HrGE',Niih_data',Niil_data'...
                ,Mh',Ml',L_data',N_data',LIC',LMIC'...
                ,UMIC',HIC',mH',mL',HnNM',HnPE',HnGE');
HumanCap = renamevars(HumanCap,["isoo","Var2","Var3","Var4","Var5","Var6"...
                               ,"Var7","Var8","Var9","Var10","Var11","Var12"...
                               ,"Var13","Var14","Var15","Var16","Var17","Var18"...
                               ,"Var19","Var20","Var21","Var22","Var23","Var24","Var25"]...
                               ,["iso","Lambda","LambdaNM","LambdaPESC","LambdaSC","Hr"...
                               ,"Hn","HrNM","HrPESC","HnSC","Niih","Niil"...
                               ,"Mh","Ml","L","N","LIC","LMIC"...
                               ,"UMIC","HIC","mH","mL","HnNM","HnPE","HnGE"]);
% Export to CSV 
writetable(HumanCap, 'matlab\output\HumanCap_Bench_25.csv', 'Delimiter',',' ,'QuoteStrings', true)


%------------------------
% Inequality results    
%------------------------ 
THEIL_Obs =[Theil;TheilA;TheilW];
THEIL_NMc =[TheilGE;TheilAGE;TheilWGE];
THEIL_NM =[TheilGEc;TheilAGEc;TheilWGEc];
THEIL_NSc =[TheilGE4;TheilAGE4;TheilWGE4];
THEIL_NS =[TheilGE4c;TheilAGE4c;TheilWGE4c];

THEIL_Bench = table(THEIL_Obs,THEIL_NMc,THEIL_NM,THEIL_NSc,THEIL_NS);
writetable(THEIL_Bench, 'matlab\output\THEIL_Bench_25.csv', 'Delimiter',',' ,'QuoteStrings', true)

toc 

