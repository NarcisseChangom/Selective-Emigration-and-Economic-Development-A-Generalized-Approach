

%-------------------------------------------------------------------------%
%               Monte carlo like calibration of z                         %
%-------------------------------------------------------------------------%
global mu zstar pol Inc_group o Lambd 
global gammahat gammahat_LIC gammahat_LMIC gammahat_UMIC gammahat_HIC gammahat_all
global zstar_lic zstar_lmic zstar_umic zstar_hic


Inc_group=zeros(1,o);
for i=1:174
    
    if LIC(i)==1
        Inc_group(i) = 1;
    end
    if LMIC(i)==1
       Inc_group(i) = 2;
    end 
    if UMIC(i)==1
        Inc_group(i) = 3;
        
    end
    if HIC(i)==1
        Inc_group(i) = 4;
    end 
    
end

% z* is the parameter governing the slope of the density function of the
% education cost/effort. z \in [0,\inf). for z=0, the distribution is
% uniform. For z>0, there are more individuals facing large education costs
% than individuals facing low education cost.

%%%%%%%%%%%%%%% SR Semi-elasticity from CDDM(2021)Estimated as part of this 
%%%%%%%%%%%%%%% paper (updated estimates) 
                gammahat=3.238;
                gammahat_LIC=3.238;
                gammahat_LMIC=3.238;
                gammahat_UMIC=0;
                gammahat_HIC=0;



%--------------------------------------

 cijh = max(0,1-(wH_data./wH_data').*((MijH_data./Niih_data).^mu));
 cijl = max(0,1-(wL_data./wL_data').*((MijL_data./Niil_data).^mu));

Wijh = (wH_data'.*(1-cijh)).^(1/mu); 
Wijl = (wL_data'.*(1-cijl)).^(1/mu); 
SumWikh = wH_data.^(1/mu) + sum((1-dums).*Wijh,1);
SumWikl = wL_data.^(1/mu) + sum((1-dums).*Wijl,1);
Lambd = SumWikh./SumWikl;
mih = sum((1-dums).*Wijh,1)./SumWikh;
mil = sum((1-dums).*Wijl,1)./SumWikl; 

%%%% COUNTERFACTUAL#1/POLICY#1 ==> cij,s when (1-cij,s)/u ==> 1-cij,sCF i.e. 
%%%% cij,sCF11 = 1-pol(1-cijs) i.e. HS and LS are affected identically
%%%% (pol is later on refers to as "beta")

%%
pol = 0.7; % ==>pol = policy. 
cijhCF11 = 1-pol*(1-cijh);
cijhCF11(cijh==0) = 0;
cijhCF11(cijhCF11>1) = 1;
cijhCF11(cijhCF11<0) = 0;
cijlCF11 = 1-pol*(1-cijl);
cijlCF11(cijl==0) = 0;
cijlCF11(cijlCF11>1) = 1;
cijlCF11(cijlCF11<0) = 0;

%%%% Calibration of Gi(z)- Gi proxies somehow the provision of public education
v=0.0:0.1:7.0;
n=length(v);
G=zeros(174,n);
gammahat_all =ones(1,n).*gammahat;
H0 = Nh_data./N_data;
h0 = Lh_data./L_data;

for j = 1:n 
    for i = 1:o
      G(i,j)= H0(i)^(1/(1+v(j)))/(1-1/Lambd(i));    
    end     
end 

%%%%% Calibration of DlnH
WijhCF11 = (wH_data'.*(1-cijhCF11)).^(1/mu); 
WijlCF11 = (wL_data'.*(1-cijlCF11)).^(1/mu); 
SumWikhCF11 = wH_data.^(1/mu) + sum((1-dums).*WijhCF11,1);
SumWiklCF11 = wL_data.^(1/mu) + sum((1-dums).*WijlCF11,1);
LambdCF11 = SumWikhCF11./SumWiklCF11;
RatioCF11 =1 - 1./LambdCF11;

for j = 1:n %generate counterfactual HC for each value of z in the range (0,oo)
 
    for i = 1:o
      Hcf11(i,j)= (G(i,j)^(1+v(j)))*RatioCF11(i)^(1+v(j));    
    end 
    
end 

mihCF11 = sum((1-dums).*WijhCF11,1)./SumWikhCF11;
milCF11 = sum((1-dums).*WijlCF11,1)./SumWiklCF11;

for j = 1:n
    for i = 1:o
      DlnH11(i,j)= log(Hcf11(i,j)) - log(H0(i)); % log deviation in HC
    end  
end 
indep11  =  (mihCF11 - milCF11)-(mih - mil); % explanotary variable 


%%%%%%%%%% COUNTERFACTUAL#2/POLICY#2. cijhCF = 1-pol(1-cijh) && cijlCF = cijl
%%
pol =0.7;
cijhCF12 = 1-pol*(1-cijh);
cijhCF12(cijh==0) = 0;
cijhCF12(cijhCF12>1) = 1;
cijhCF12(cijhCF12<0) = 0;
cijlCF12 = cijl;

%%%%% Calibration of DlnH
WijhCF12 = (wH_data'.*(1-cijhCF12)).^(1/mu); 
WijlCF12 = (wL_data'.*(1-cijlCF12)).^(1/mu); 
SumWikhCF12 = wH_data.^(1/mu) + sum((1-dums).*WijhCF12,1);
SumWiklCF12 = wL_data.^(1/mu) + sum((1-dums).*WijlCF12,1);
LambdCF12 = SumWikhCF12./SumWiklCF12;
RatioCF12 =1 - 1./LambdCF12;

for j = 1:n %generate counterfactual HC for each value of z in the range (0,oo)
 
    for i = 1:o
      Hcf12(i,j)= (G(i,j)^(1+v(j)))*RatioCF12(i)^(1+v(j));    
    end 
    
end 

mihCF12 = sum((1-dums).*WijhCF12,1)./SumWikhCF12;
milCF12 = sum((1-dums).*WijlCF12,1)./SumWiklCF12; 
for j = 1:n
    for i = 1:o
      DlnH12(i,j)= log(Hcf12(i,j)) - log(H0(i)); % log deviation in HC
    end  
end 
indep12  =  (mihCF12 - milCF12)-(mih - mil); % explanotary variable 
  

%%%%%%%%%% COUNTERFACTUAL#3/POLICY#3 - cijlCF = 1-pol(1-cijl) & cijhCF = cijh
%%
pol =0.7;
cijhCF13 = cijh;
cijlCF13 = 1-pol*(1-cijl);
cijlCF13(cijl==0) = 0;
cijlCF13(cijlCF13>1) = 1;
cijlCF13(cijlCF13<0) = 0;

%%%%% Calibration of DlnH
WijhCF13 = (wH_data'.*(1-cijhCF13)).^(1/mu); 
WijlCF13 = (wL_data'.*(1-cijlCF13)).^(1/mu); 
SumWikhCF13 = wH_data.^(1/mu) + sum((1-dums).*WijhCF13,1);
SumWiklCF13 = wL_data.^(1/mu) + sum((1-dums).*WijlCF13,1);
LambdCF13 = SumWikhCF13./SumWiklCF13;
RatioCF13 =1 - 1./LambdCF13;

for j = 1:n %generate counterfactual HC for each value of z in the range (0,oo)
 
    for i = 1:o
      Hcf13(i,j)= (G(i,j)^(1+v(j)))*RatioCF13(i)^(1+v(j));    
    end 
    
end 

mihCF13 = sum((1-dums).*WijhCF13,1)./SumWikhCF13;
milCF13 = sum((1-dums).*WijlCF13,1)./SumWiklCF13; 
for j = 1:n
    for i = 1:o
      DlnH13(i,j)= log(Hcf13(i,j)) - log(H0(i)); % log deviation in HC
    end  
end
indep13  =  (mihCF13 - milCF13)-(mih - mil); % explanotary variable 


%%%%%%%%%% COUNTERFACTUAL #4\POLICY#4 - u(1-cijs)=> i.e. 1.3(1-cijs)
%%
pol = 1.3;
cijhCF14 = 1-pol*(1-cijh);
cijhCF14(cijh==0) = 0;
cijhCF14(cijhCF14>1) = 1;
cijhCF14(cijhCF14<0) = 0;
cijlCF14 = 1-pol*(1-cijl);
cijlCF14(cijl==0) = 0;
cijlCF14(cijlCF14>1) = 1;
cijlCF14(cijlCF14<0) = 0;

%%%%% Calibration of DlnH
WijhCF14 = (wH_data'.*(1-cijhCF14)).^(1/mu); 
WijlCF14 = (wL_data'.*(1-cijlCF14)).^(1/mu); 
SumWikhCF14 = wH_data.^(1/mu) + sum((1-dums).*WijhCF14,1);
SumWiklCF14 = wL_data.^(1/mu) + sum((1-dums).*WijlCF14,1);
LambdCF14 = SumWikhCF14./SumWiklCF14;
RatioCF14 =1 - 1./LambdCF14;

for j = 1:n %generate counterfactual HC for each value of z in the range (0,oo)
 
    for i = 1:o
      Hcf14(i,j)= (G(i,j)^(1+v(j)))*RatioCF14(i)^(1+v(j));    
    end 
    
end 

mihCF14 = sum((1-dums).*WijhCF14,1)./SumWikhCF14;
milCF14 = sum((1-dums).*WijlCF14,1)./SumWiklCF14; 
for j = 1:n
    for i = 1:o
      DlnH14(i,j)= log(Hcf14(i,j)) - log(H0(i)); % log deviation in HC
    end  
end
indep14  =  (mihCF14 - milCF14)-(mih - mil); % explanotary variable 


%%%%%%%%%% COUNTERFACTUAL #5/POLICY#5 => u(1-cijs): i.e. 1.3(1-cijh)
pol = 1.3; %  ==>pol = policy
cijhCF15 = 1-pol*(1-cijh);
cijhCF15(cijh==0) = 0;
cijhCF15(cijhCF15>1) = 1;
cijhCF15(cijhCF15<0) = 0;
cijlCF15 = cijl;

%%%%% Calibration of DlnH
WijhCF15 = (wH_data'.*(1-cijhCF15)).^(1/mu); 
WijlCF15 = (wL_data'.*(1-cijlCF15)).^(1/mu); 
SumWikhCF15 = wH_data.^(1/mu) + sum((1-dums).*WijhCF15,1);
SumWiklCF15 = wL_data.^(1/mu) + sum((1-dums).*WijlCF15,1);
LambdCF15 = SumWikhCF15./SumWiklCF15;
RatioCF15 =1 - 1./LambdCF15;

for j = 1:n %generate counterfactual HC for each value of z in the range (0,oo)
 
    for i = 1:o
      Hcf15(i,j)= (G(i,j)^(1+v(j)))*RatioCF15(i)^(1+v(j));    
    end 
    
end 

mihCF15 = sum((1-dums).*WijhCF15,1)./SumWikhCF15;
milCF15 = sum((1-dums).*WijlCF15,1)./SumWiklCF15; 
for j = 1:n
    for i = 1:o
      DlnH15(i,j)= log(Hcf15(i,j)) - log(H0(i)); % log deviation in HC
    end  
end
indep15  =  (mihCF15 - milCF15)-(mih - mil); % explanotary variable 


%%%%%%%%%% COUNTERFACTUAL #6 /POLICY #6 u(1-cijs) => 1.3(1-cijl)
%%
pol = 1.3; % i.e.  ==>pol = policy
cijhCF16 = cijh;
cijlCF16 = 1-pol*(1-cijl);
cijlCF16(cijl==0) = 0;
cijlCF16(cijlCF16>1) = 1;
cijlCF16(cijlCF16<0) = 0;


%%%%% Calibration of DlnH
WijhCF16 = (wH_data'.*(1-cijhCF16)).^(1/mu); 
WijlCF16 = (wL_data'.*(1-cijlCF16)).^(1/mu); 
SumWikhCF16 = wH_data.^(1/mu) + sum((1-dums).*WijhCF16,1);
SumWiklCF16 = wL_data.^(1/mu) + sum((1-dums).*WijlCF16,1);
LambdCF16 = SumWikhCF16./SumWiklCF16;
RatioCF16 =1 - 1./LambdCF16;

for j = 1:n %generate counterfactual HC for each value of z in the range (0,oo)
 
    for i = 1:o
      Hcf16(i,j)= (G(i,j)^(1+v(j)))*RatioCF16(i)^(1+v(j));    
    end 
    
end 

mihCF16 = sum((1-dums).*WijhCF16,1)./SumWikhCF16;
milCF16 = sum((1-dums).*WijlCF16,1)./SumWiklCF16; 
for j = 1:n
    for i = 1:o
      DlnH16(i,j)= log(Hcf16(i,j)) - log(H0(i)); % log deviation in HC
    end  
end
indep16  =  (mihCF16 - milCF16)-(mih - mil); % explanotary variable 


%%%%%%%%%% COUNTERFACTUAL #7 /POLICY #6 u(1-cijs) => 0(1-cijs)
cijhCF17 = ones(o,d);
cijlCF17 = ones(o,d);


%%%%% Calibration of DlnH
WijhCF17 = (wH_data'.*(1-cijhCF17)).^(1/mu); 
WijlCF17 = (wL_data'.*(1-cijlCF17)).^(1/mu); 
SumWikhCF17 = wH_data.^(1/mu) + sum((1-dums).*WijhCF17,1);
SumWiklCF17 = wL_data.^(1/mu) + sum((1-dums).*WijlCF17,1);
LambdCF17 = SumWikhCF17./SumWiklCF17;
RatioCF17 =1 - 1./LambdCF17;

for j = 1:n %generate counterfactual HC for each value of z in the range (0,oo)
 
    for i = 1:o
      Hcf17(i,j)= (G(i,j)^(1+v(j)))*RatioCF17(i)^(1+v(j));    
    end 
    
end 

mihCF17 = sum((1-dums).*WijhCF17,1)./SumWikhCF17;
milCF17 = sum((1-dums).*WijlCF17,1)./SumWiklCF17; 
for j = 1:n
    for i = 1:o
      DlnH17(i,j)= log(Hcf17(i,j)) - log(H0(i)); % log deviation in HC
    end  
end
indep17  =  (mihCF17 - milCF17)-(mih - mil); % explanotary variable 


%%------------------------  ESTIMATING z --------------------------------%%
%%-----------------------------------------------------------------------%%

%%%%% country specific theoretical semi-elasticities
for j=1:n    
    for i=1:o
      eT1(i,j)=DlnH11(i,j)/indep11(i); % eT is the theoretical elasticity for each country in the sample  
      eT2(i,j)=DlnH12(i,j)/indep12(i);
      eT3(i,j)=DlnH13(i,j)/indep13(i);
      eT4(i,j)=DlnH14(i,j)/indep14(i);
      eT5(i,j)=DlnH15(i,j)/indep15(i);
      eT6(i,j)=DlnH16(i,j)/indep16(i);
      eT7(i,j)=DlnH17(i,j)/indep17(i);
    end    
end

%%%%% benchmark values of elasticities (literature)
Ei_hat =zeros(1,174);
for i=1:174
    if LIC(i)==1
        Ei_hat(i) = gammahat_LIC;
    end
    
     if LMIC(i)==1
         Ei_hat(i) = gammahat_LMIC;
     end
         if UMIC(i)==1
             Ei_hat(i) = gammahat_UMIC;
         end 
            if HIC(i)==1
                 Ei_hat(i) = gammahat_HIC;
            end 
end
 
%%%%%%%% Gap between theoritical country specific elasticities and
%%%%%%%% benchmark values and RSS
for j=1:n   
    for i=1:o
    eTminusEi_hat1(i,j) = (eT1(i,j) - Ei_hat(i))*(1-HIC(i));   
    eTminusEi_hat2(i,j) = (eT2(i,j) - Ei_hat(i))*(1-HIC(i));
    eTminusEi_hat3(i,j) = (eT3(i,j) - Ei_hat(i))*(1-HIC(i));
    eTminusEi_hat4(i,j) = (eT4(i,j) - Ei_hat(i))*(1-HIC(i));
    eTminusEi_hat5(i,j) = (eT5(i,j) - Ei_hat(i))*(1-HIC(i));
    eTminusEi_hat6(i,j) = (eT6(i,j) - Ei_hat(i))*(1-HIC(i));
    eTminusEi_hat7(i,j) = (eT7(i,j) - Ei_hat(i))*(1-HIC(i));
    end 
 MSE1(j)= (sum(eTminusEi_hat1(:,j))^2)/(sum(LIC(1,:))+sum(LMIC(1,:))+sum(UMIC(1,:))); 
 MSE2(j)= (sum(eTminusEi_hat2(:,j))^2)/(sum(LIC(1,:))+sum(LMIC(1,:))+sum(UMIC(1,:)));
 MSE3(j)= (sum(eTminusEi_hat3(:,j))^2)/(sum(LIC(1,:))+sum(LMIC(1,:))+sum(UMIC(1,:)));
 MSE4(j)= (sum(eTminusEi_hat4(:,j))^2)/(sum(LIC(1,:))+sum(LMIC(1,:))+sum(UMIC(1,:)));
 MSE5(j)= (sum(eTminusEi_hat5(:,j))^2)/(sum(LIC(1,:))+sum(LMIC(1,:))+sum(UMIC(1,:))); 
 MSE6(j)= (sum(eTminusEi_hat6(:,j))^2)/(sum(LIC(1,:))+sum(LMIC(1,:))+sum(UMIC(1,:)));
 MSE7(j)= (sum(eTminusEi_hat7(:,j))^2)/(sum(LIC(1,:))+sum(LMIC(1,:))+sum(UMIC(1,:))); 
end

%%%%%%%% BY INCOME GROUP

%%%%% 1. LOW INCOME COUNTRIES
for j=1:n   
    for i=1:o
    eTminusEi_hat1(i,j) = (eT1(i,j) - Ei_hat(i))*LIC(i);   
    eTminusEi_hat2(i,j) = (eT2(i,j) - Ei_hat(i))*LIC(i);
    eTminusEi_hat3(i,j) = (eT3(i,j) - Ei_hat(i))*LIC(i);
    eTminusEi_hat4(i,j) = (eT4(i,j) - Ei_hat(i))*LIC(i);
    eTminusEi_hat5(i,j) = (eT5(i,j) - Ei_hat(i))*LIC(i);
    eTminusEi_hat6(i,j) = (eT6(i,j) - Ei_hat(i))*LIC(i);
    eTminusEi_hat7(i,j) = (eT7(i,j) - Ei_hat(i))*LIC(i);
    end 
 MSE1LIC(j)= (sum(eTminusEi_hat1(:,j))^2)/sum(LIC(1,:)); 
 MSE2LIC(j)= (sum(eTminusEi_hat2(:,j))^2)/sum(LIC(1,:));
 MSE3LIC(j)= (sum(eTminusEi_hat3(:,j))^2)/sum(LIC(1,:));
 MSE4LIC(j)= (sum(eTminusEi_hat4(:,j))^2)/sum(LIC(1,:));
 MSE5LIC(j)= (sum(eTminusEi_hat5(:,j))^2)/sum(LIC(1,:)); 
 MSE6LIC(j)= (sum(eTminusEi_hat6(:,j))^2)/sum(LIC(1,:));  
 MSE7LIC(j)= (sum(eTminusEi_hat7(:,j))^2)/sum(LIC(1,:));  
end

%%%%% 2. LOWER MIDDLE INCOME COUNTRIES

for j=1:n
    
    for i=1:174
    eTminusEi_hat1lmic(i,j) = (eT1(i,j) - Ei_hat(i))*LMIC(i);   
    eTminusEi_hat2lmic(i,j) = (eT2(i,j) - Ei_hat(i))*LMIC(i);
    eTminusEi_hat3lmic(i,j) = (eT3(i,j) - Ei_hat(i))*LMIC(i);
    eTminusEi_hat4lmic(i,j) = (eT4(i,j) - Ei_hat(i))*LMIC(i);
    eTminusEi_hat5lmic(i,j) = (eT5(i,j) - Ei_hat(i))*LMIC(i);
    eTminusEi_hat6lmic(i,j) = (eT6(i,j) - Ei_hat(i))*LMIC(i);
     eTminusEi_hat7lmic(i,j) = (eT7(i,j) - Ei_hat(i))*LMIC(i);
    end 
 MSE1LMIC(j)= (sum(eTminusEi_hat1lmic(:,j))^2)/sum(LMIC(1,:)); 
 MSE2LMIC(j)= (sum(eTminusEi_hat2lmic(:,j))^2)/sum(LMIC(1,:));
 MSE3LMIC(j)= (sum(eTminusEi_hat3lmic(:,j))^2)/sum(LMIC(1,:));
 MSE4LMIC(j)= (sum(eTminusEi_hat4lmic(:,j))^2)/sum(LMIC(1,:));
 MSE5LMIC(j)= (sum(eTminusEi_hat5lmic(:,j))^2)/sum(LMIC(1,:)); 
 MSE6LMIC(j)= (sum(eTminusEi_hat6lmic(:,j))^2)/sum(LMIC(1,:)); 
 MSE7LMIC(j)= (sum(eTminusEi_hat7lmic(:,j))^2)/sum(LMIC(1,:));
end

%%%%%% 3. UPPER MIDDLE INCOME COUNTRIES


for j=1:n
    
    for i=1:174
    eTminusEi_hat1umic(i,j) = (eT1(i,j) - Ei_hat(i))*UMIC(i);   
    eTminusEi_hat2umic(i,j) = (eT2(i,j) - Ei_hat(i))*UMIC(i);
    eTminusEi_hat3umic(i,j) = (eT3(i,j) - Ei_hat(i))*UMIC(i);
    eTminusEi_hat4umic(i,j) = (eT4(i,j) - Ei_hat(i))*UMIC(i);
    eTminusEi_hat5umic(i,j) = (eT5(i,j) - Ei_hat(i))*UMIC(i);
    eTminusEi_hat6umic(i,j) = (eT6(i,j) - Ei_hat(i))*UMIC(i);
    eTminusEi_hat7umic(i,j) = (eT7(i,j) - Ei_hat(i))*UMIC(i);
    end 
 MSE1UMIC(j)= (sum(eTminusEi_hat1umic(:,j))^2)/sum(UMIC(1,:)); 
 MSE2UMIC(j)= (sum(eTminusEi_hat2umic(:,j))^2)/sum(UMIC(1,:));
 MSE3UMIC(j)= (sum(eTminusEi_hat3umic(:,j))^2)/sum(UMIC(1,:));
 MSE4UMIC(j)= (sum(eTminusEi_hat4umic(:,j))^2)/sum(UMIC(1,:));
 MSE5UMIC(j)= (sum(eTminusEi_hat5umic(:,j))^2)/sum(UMIC(1,:)); 
 MSE6UMIC(j)= (sum(eTminusEi_hat6umic(:,j))^2)/sum(UMIC(1,:));  
 MSE7UMIC(j)= (sum(eTminusEi_hat7umic(:,j))^2)/sum(UMIC(1,:));  
end



%%%%%% 4. HIGH INCOME COUNTRIES


for j=1:n
    
    for i=1:174
    eTminusEi_hat1hic(i,j) = (eT1(i,j) - Ei_hat(i))*HIC(i);   
    eTminusEi_hat2hic(i,j) = (eT2(i,j) - Ei_hat(i))*HIC(i);
    eTminusEi_hat3hic(i,j) = (eT3(i,j) - Ei_hat(i))*HIC(i);
    eTminusEi_hat4hic(i,j) = (eT4(i,j) - Ei_hat(i))*HIC(i);
    eTminusEi_hat5hic(i,j) = (eT5(i,j) - Ei_hat(i))*HIC(i);
    eTminusEi_hat6hic(i,j) = (eT6(i,j) - Ei_hat(i))*HIC(i);
eTminusEi_hat7hic(i,j) = (eT7(i,j) - Ei_hat(i))*HIC(i);
    end 
 MSE1HIC(j)= (sum(eTminusEi_hat1hic(:,j))^2)/sum(HIC(1,:)); 
 MSE2HIC(j)= (sum(eTminusEi_hat2hic(:,j))^2)/sum(HIC(1,:));
 MSE3HIC(j)= (sum(eTminusEi_hat3hic(:,j))^2)/sum(HIC(1,:));
 MSE4HIC(j)= (sum(eTminusEi_hat4hic(:,j))^2)/sum(HIC(1,:));
 MSE5HIC(j)= (sum(eTminusEi_hat5hic(:,j))^2)/sum(HIC(1,:)); 
 MSE6HIC(j)= (sum(eTminusEi_hat6hic(:,j))^2)/sum(HIC(1,:));  
  MSE7HIC(j)= (sum(eTminusEi_hat7hic(:,j))^2)/sum(HIC(1,:)); 
end


%%%%%%%%%%%%%%%% Combined graph for the optimal z by income group %%%%%%%%%

tiledlayout(2,2)  
nexttile 
plots=plot(v,MSE7LIC,'-- k','LineWidth',1.3);
xticks(0:1:8);
names={'1-c_{ij,s} =0'};
xlabel('z') 
ylabel('RSS (z)') 


nexttile 
plots=plot(v,MSE7LMIC,'-- k','LineWidth',1.3);
xticks(0:1:8);
names={'1-c_{ij,s} =0'};
xlabel('z') 
ylabel('RSS (z)') 


nexttile 
plots=plot(v,MSE7UMIC,'-- k','LineWidth',1.3);
xticks(0:1:8);
names={'1-c_{ij,s} =0'};
xlabel('z') 
ylabel('RSS (z)') 

nexttile 
plots=plot(v,MSE7HIC,'-- k','LineWidth',1.3);
xticks(0:1:8);
names={'1-c_{ij,s} =0'};
xlabel('z') 
ylabel('RSS (z)') 

%saveas(gcf,'figures\Figure_B1.png')
saveas(gcf, '..\main\figures\Figure_B1.png');    
