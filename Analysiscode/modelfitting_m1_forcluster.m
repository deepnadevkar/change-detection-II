function modelfitting_m1(sind)
tic
load datasets_m1
sessions = datasets;

conddata = sessions{sind}(:,1);
deltadata = sessions{sind}(:,2);
corrdata = sessions{sind}(:,3);

condvec = 1:4;
deltavec = 10:10:90;

for condind = 1:length(condvec)
    cond = condvec(condind);
    idx = find(conddata == cond);
    corrdataidx_1_cond(condind)= sum(corrdata(idx) ==1);
    corrdataidx_0_cond(condind)= sum(corrdata(idx) ==0);
    perfdata1(condind) = mean(corrdata(conddata == cond)); 
    for deltaind = 1:length(deltavec)
        delta = deltavec(deltaind);
        idx2 = find(conddata == cond & deltadata == delta);
        corrdataidx_1(condind, deltaind) = sum(corrdata(idx2)==1);
        corrdataidx_0(condind, deltaind)= sum(corrdata(idx2)==0);
        perfdata2(condind,deltaind) = mean(corrdata(conddata == cond & deltadata == delta)); 
    end
end

%% Fitting models (by maximizing likelihood)   

%VSO model fit
load modelpred_VSO

loglike_VSO = bsxfun(@times, corrdataidx_1, log(perfmodel_VSO))+ bsxfun(@times, corrdataidx_0, log(1 - perfmodel_VSO));
totalloglike_VSO = squeeze(sum(sum(loglike_VSO,2),1)); 
[LL_VSO, idx2] = max(totalloglike_VSO(:));
[Jbarhmax_idx_VSO, Jbarlmax_idx_VSO,taumax_idx_VSO] = ind2sub(size(totalloglike_VSO),idx2);
pars_VSO = [Jbarhvec(Jbarhmax_idx_VSO) Jbarlvec(Jbarlmax_idx_VSO) tauvec(taumax_idx_VSO)];
perfmodel1_VSO = squeeze(mean(perfmodel_VSO,2));
VSO1 = squeeze(perfmodel1_VSO(:,Jbarhmax_idx_VSO,Jbarlmax_idx_VSO,taumax_idx_VSO));
VSO2 = squeeze(perfmodel_VSO(:,:,Jbarhmax_idx_VSO,Jbarlmax_idx_VSO,taumax_idx_VSO));

Jbarhmax = max(Jbarhvec); Jbarhsteps = Jbarhvec(1,2); Jbarlmax = max(Jbarlvec); Jbarlsteps = Jbarlvec(1,2); taumax= max(tauvec); tausteps = tauvec(1,2)-tauvec(1,1);  
BMC_VSO = log((Jbarhsteps/Jbarhmax)*(Jbarlsteps/Jbarlmax)*(tausteps/taumax))+LL_VSO+log(sum(exp(totalloglike_VSO(:)-LL_VSO)));
BIC_VSO = ((-2*LL_VSO + 3*log(length(conddata))))/-2; 
AIC_VSO = ((-2*LL_VSO + 3*2))/-2; 
AICc_VSO = AIC_VSO +(24/(length(conddata))-2);


RMSE_SS_VSO = sqrt(mean((perfdata1(:)- VSO1(:)).^2));
RMSE_CM_VSO = sqrt(mean((perfdata2(:)- VSO2(:)).^2));
Rsquared_SS_VSO = 1-(sum((perfdata1(:)- VSO1(:)).^2))/sum((perfdata1(:)- mean(perfdata1(:))).^2);
Rsquared_CM_VSO = 1-(sum((perfdata2(:)- VSO2(:)).^2))/sum((perfdata2(:)- mean(perfdata2(:))).^2);
Rs_VSO = [RMSE_SS_VSO RMSE_CM_VSO Rsquared_SS_VSO Rsquared_CM_VSO];

fprintf('Done VSO model\n');
clear perfmodel_VSO;

%VEO model fit
load modelpred_VEO

loglike_VEO = bsxfun(@times, corrdataidx_1, log(perfmodel_VEO))+ bsxfun(@times, corrdataidx_0, log(1 - perfmodel_VEO));
totalloglike_VEO = squeeze(sum(sum(loglike_VEO,2),1)); 
[LL_VEO, idx2] = max(totalloglike_VEO(:));
[Jbarhmax_idx_VEO, Jbarlmax_idx_VEO,taumax_idx_VEO] = ind2sub(size(totalloglike_VEO),idx2);
pars_VEO = [Jbarhvec(Jbarhmax_idx_VEO) Jbarlvec(Jbarlmax_idx_VEO) tauvec(taumax_idx_VEO)];
perfmodel1_VEO = squeeze(mean(perfmodel_VEO,2));
VEO1 = squeeze(perfmodel1_VEO(:,Jbarhmax_idx_VEO,Jbarlmax_idx_VEO,taumax_idx_VEO));
VEO2 = squeeze(perfmodel_VEO(:,:,Jbarhmax_idx_VEO,Jbarlmax_idx_VEO,taumax_idx_VEO));

Jbarhmax = max(Jbarhvec); Jbarhsteps = Jbarhvec(1,2); Jbarlmax = max(Jbarlvec); Jbarlsteps = Jbarlvec(1,2); taumax= max(tauvec); tausteps = tauvec(1,2)-tauvec(1,1);  
BMC_VEO = log((Jbarhsteps/Jbarhmax)*(Jbarlsteps/Jbarlmax)*(tausteps/taumax))+LL_VEO+log(sum(exp(totalloglike_VEO(:)-LL_VEO)));
BIC_VEO = ((-2*LL_VEO + 3*log(length(conddata))))/-2; 
AIC_VEO = ((-2*LL_VEO + 3*2))/-2; 
AICc_VEO = AIC_VEO +(24/(length(conddata))-2);

RMSE_SS_VEO = sqrt(mean((perfdata1(:)- VEO1(:)).^2));
RMSE_CM_VEO = sqrt(mean((perfdata2(:)- VEO2(:)).^2));
Rsquared_SS_VEO = 1-(sum((perfdata1(:)- VEO1(:)).^2))/sum((perfdata1(:)- mean(perfdata1(:))).^2);
Rsquared_CM_VEO = 1-(sum((perfdata2(:)- VEO2(:)).^2))/sum((perfdata2(:)- mean(perfdata2(:))).^2);
Rs_VEO = [RMSE_SS_VEO RMSE_CM_VEO Rsquared_SS_VEO Rsquared_CM_VEO];

fprintf('Done VEO model\n');
clear perfmodel_VEO;

%VVO model fit
load modelpred_VVO

loglike_VVO = bsxfun(@times, corrdataidx_1, log(perfmodel_VVO))+ bsxfun(@times, corrdataidx_0, log(1 - perfmodel_VVO));
totalloglike_VVO = squeeze(sum(sum(loglike_VVO,2),1)); 
[LL_VVO, idx2] = max(totalloglike_VVO(:));
[Jbarhmax_idx_VVO, Jbarlmax_idx_VVO,taumax_idx_VVO] = ind2sub(size(totalloglike_VVO),idx2);
pars_VVO = [Jbarhvec(Jbarhmax_idx_VVO) Jbarlvec(Jbarlmax_idx_VVO) tauvec(taumax_idx_VVO)];
perfmodel1_VVO = squeeze(mean(perfmodel_VVO,2));
VVO1 = squeeze(perfmodel1_VVO(:,Jbarhmax_idx_VVO,Jbarlmax_idx_VVO,taumax_idx_VVO));
VVO2 = squeeze(perfmodel_VVO(:,:,Jbarhmax_idx_VVO,Jbarlmax_idx_VVO,taumax_idx_VVO));

Jbarhmax = max(Jbarhvec); Jbarhsteps = Jbarhvec(1,2); Jbarlmax = max(Jbarlvec); Jbarlsteps = Jbarlvec(1,2); taumax= max(tauvec); tausteps = tauvec(1,2)-tauvec(1,1);  
BMC_VVO = log((Jbarhsteps/Jbarhmax)*(Jbarlsteps/Jbarlmax)*(tausteps/taumax))+LL_VVO+log(sum(exp(totalloglike_VVO(:)-LL_VVO)));
BIC_VVO = ((-2*LL_VVO + 3*log(length(conddata))))/-2; 
AIC_VVO = ((-2*LL_VVO + 3*2))/-2; 
AICc_VVO = AIC_VVO +(24/(length(conddata))-2);

RMSE_SS_VVO = sqrt(mean((perfdata1(:)- VVO1(:)).^2));
RMSE_CM_VVO = sqrt(mean((perfdata2(:)- VVO2(:)).^2));
Rsquared_SS_VVO = 1-(sum((perfdata1(:)- VVO1(:)).^2))/sum((perfdata1(:)- mean(perfdata1(:))).^2);
Rsquared_CM_VVO = 1-(sum((perfdata2(:)- VVO2(:)).^2))/sum((perfdata2(:)- mean(perfdata2(:))).^2);
Rs_VVO = [RMSE_SS_VVO RMSE_CM_VVO Rsquared_SS_VVO Rsquared_CM_VVO];

fprintf('Done VVO model\n');
clear perfmodel_VVO;

%Report all
LL = [LL_VSO LL_VEO LL_VVO];
BMC = [BMC_VSO BMC_VEO BMC_VVO];
BIC = [BIC_VSO BIC_VEO BIC_VVO];
AIC = [AIC_VSO AIC_VEO AIC_VVO];
AICc = [AICc_VSO AICc_VEO AICc_VVO];


filename = strcat('modelfitsM1_',num2str(sind))
save(filename, 'VSO1', 'VSO2', 'VEO1', 'VEO2', 'VVO1', 'VVO2', 'pars_VSO', 'pars_VEO', 'pars_VVO', 'LL', 'BMC', 'BIC', 'AIC', 'AICc', 'Rs_VSO', 'Rs_VEO', 'Rs_VVO');
toc


