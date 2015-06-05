clear all close all;

for ii = 1:100
    fname = sprintf('/Local/home/deepna/Old Lab Backup/Wright Lab/Desk Computer/CDTIICode - 4.4.2015 - withlapse/modelfits_m1/modelfits_m1_%i.mat', ii);
    load(fname);
    AIC_all(ii,:) = AIC;
    AICc_all(ii,:) = AICc;
    BIC_all(ii,:) = BIC;
    BMC_all(ii,:) = BMC;
    LL_all(ii,:) = LL;
    VSO1_all(ii,:) = VSO1;
    VSO2_all(ii,:,:) = VSO2;
    VEO1_all(ii,:) = VEO1;
    VEO2_all(ii,:,:) = VEO2;
    VVO1_all(ii,:)= VVO1;
    VVO2_all(ii,:,:) = VVO2;
    pars_VEO_all(ii,:) = pars_VEO;
    pars_VSO_all(ii,:) = pars_VSO;
    pars_VVO_all(ii,:) = pars_VVO;
    Rs_VEO_all(ii,:) = Rs_VEO;
    Rs_VSO_all(ii,:) = Rs_VSO;
    Rs_VVO_all(ii,:) = Rs_VVO;
end

save modelfitting_m1