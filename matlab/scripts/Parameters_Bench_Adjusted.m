global mu sigma psi kappa epsilon eta lambda1 rho mbar o zstar_lic zstar_lmic zstar zstar_umic zstar_hic

%% Parameters 
%mu = 1/0.7; 
mu=0.7; sigma =2; psi = (sigma - 1)/sigma; Pi=1;  kappa = 0.1; epsilon =0.1; 
eta = 0.056;  lambda1 = 8;  rho = 0.032; mbar= 0.1; 
zstar_lic =23.4; zstar_lmic=15.9; zstar_umic=0.0; zstar_hic=0.0;

% Define zstar across income groups  
LIC =zeros(1,o); LMIC=zeros(1,o); UMIC=zeros(1,o); HIC=zeros(1,o); zstar=zeros(1,174);

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



for i=1:o
    if LIC(i)==1
        zstar(i)=zstar_lic;
    end
    if LMIC(i)==1
        zstar(i)=zstar_lmic;
    end
    if UMIC(i)==1
        zstar(i)=zstar_umic;
    end 
    if HIC(i)==1
        zstar(i)=zstar_hic;
    end
end
