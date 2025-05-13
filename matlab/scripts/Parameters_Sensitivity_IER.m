global mu sigma psi kappa epsilon eta lambda1 rho mbar o zstar_lic zstar_lmic zstar zstar_umic zstar_hic

%% Parameters 
mu = 1/0.7; sigma =2; psi = (sigma - 1)/sigma; Pi=1;  kappa = 0.1; epsilon =0.1; 
eta = 0.056;  lambda1 = 8;  rho = 0.032; mbar= 0.1; 
zstar_lic =5.3; zstar_lmic=3.8; zstar_umic=0.0; zstar_hic=0.0;

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
