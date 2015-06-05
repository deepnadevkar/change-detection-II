clear all; close all;

load datasets_m2
sessions = datasets;

load modelfitting_m2 
    
condvec = 1:4;
deltavec = 10:10:90;

perfdata1 = zeros(length(sessions),length(condvec));
perfdata2  = zeros(length(sessions),length(condvec),length(deltavec));

%% Loading data and plotting performance across sessions

for sind = 1:1:length(sessions)
    sind
    
    conddata = sessions{sind}(:,1);
    deltadata = sessions{sind}(:,2);
    corrdata = sessions{sind}(:,3);
    
    for condind = 1:length(condvec)
        cond = condvec(condind);
        perfdata1(sind,condind) = mean(corrdata(conddata == cond));
        for deltaind = 1:length(deltavec)
            delta = deltavec(deltaind);
            perfdata2(sind,condind,deltaind) = mean(corrdata(conddata == cond & deltadata == delta));
        end
    end
end 
    
% %criterion diffs from VVO model
% AIC_diff = [AIC_all(:,1)-AIC_all(:,3) AIC_all(:,2)-AIC_all(:,3)];
% AICc_diff = [AICc_all(:,1)-AICc_all(:,3) AICc_all(:,2)-AICc_all(:,3)];
% BIC_diff = [BIC_all(:,1)-BIC_all(:,3) BIC_all(:,2)-BIC_all(:,3)];
% BMC_diff = [BMC_all(:,1)-BMC_all(:,3) BMC_all(:,2)-BMC_all(:,3)];

%criterion diffs from VEO model
AIC_diff = [AIC_all(:,1)-AIC_all(:,2) AIC_all(:,3)-AIC_all(:,2)];
AICc_diff = [AICc_all(:,1)-AICc_all(:,2) AICc_all(:,3)-AICc_all(:,2)];
BIC_diff = [BIC_all(:,1)-BIC_all(:,2) BIC_all(:,3)-BIC_all(:,2)];
BMC_diff = [BMC_all(:,1)-BMC_all(:,2) BMC_all(:,3)-BMC_all(:,2)];

%% means
datacond_mean = nanmean(perfdata1,1);
datadelta_mean = squeeze(nanmean(perfdata2,1));
VSO1_mean = mean(VSO1_all,1);
VSO2_mean = squeeze(mean(VSO2_all,1));
VEO1_mean = mean(VEO1_all,1);
VEO2_mean = squeeze(mean(VEO2_all,1));
VVO1_mean = mean(VVO1_all,1);
VVO2_mean = squeeze(mean(VVO2_all,1));
pars_VSO_mean = mean(pars_VSO_all,1);
pars_VEO_mean = mean(pars_VEO_all,1);
pars_VVO_mean = mean(pars_VVO_all,1);
BMC_mean = mean(BMC_all,1);
BIC_mean = mean(BIC_all,1);
AIC_mean = mean(AIC_all,1);
AICc_mean = mean(AICc_all,1);
Rs_VSO_mean = mean(Rs_VSO_all,1);
Rs_VEO_mean = mean(Rs_VEO_all,1);
Rs_VVO_mean = mean(Rs_VVO_all,1);
AIC_meandiff = mean(AIC_diff,1);
AICc_meandiff = mean(AICc_diff,1);
BIC_meandiff = mean(BIC_diff,1);
BMC_meandiff = mean(BMC_diff,1);


%% standard errors

datacond_ste = nanstd(perfdata1,0,1);
datadelta_ste = squeeze(nanstd(perfdata2,0,1));
VSO1_ste = std(VSO1_all,0,1);
VSO2_ste = squeeze(std(VSO2_all,0,1));
VEO1_ste = std(VEO1_all,0,1);
VEO2_ste = squeeze(std(VEO2_all,0,1));
VVO1_ste = std(VVO1_all,0,1);
VVO2_ste = squeeze(std(VVO2_all,0,1));
pars_VSO_ste = std(pars_VSO_all,0,1);
pars_VEO_ste = std(pars_VEO_all,0,1);
pars_VVO_ste = std(pars_VVO_all,0,1);
BMC_ste = std(BMC_all,0,1);
BIC_ste = std(BIC_all,0,1);
AIC_ste = std(AIC_all,0,1);
AICc_ste = std(AICc_all,0,1);
BMC_diffste = std(BMC_diff,0,1);
BIC_diffste = std(BIC_diff,0,1);
AIC_diffste = std(AIC_diff,0,1);
AICc_diffste = std(AICc_diff,0,1);
Rs_VSO_ste = std(Rs_VSO_all,0,1);
Rs_VEO_ste = std(Rs_VEO_all,0,1);
Rs_VVO_ste = std(Rs_VVO_all,0,1);

% behavioral plots
co = [0 0 1;
      0 0.5 0;
      1 0 0;
      0 0.75 0.75;
      0.75 0 0.75;
      0.75 0.75 0;
      0.25 0.25 0.25];

figure;
set(gca,'FontSize',11,'FontName','Arial');
xfac = 1;
yfac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac yfac]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac yfac]);
hold on;
set(gca,'YTick',0.2:0.2:1.0)
set(gca,'XTick',0:30:90)
errorbar(repmat(deltavec,4,1)', datadelta_mean', datadelta_ste','-o'); hold on;
xlabel('Change magnitude'); ylabel('Proportion correct');
axis([7 93 0.2 1.0]);
legend(strcat('Condition= ',int2str(condvec')), 4, 'Location', 'Best');

%% Condition plots

%VSO
figure; %graphs performance across conditions with error bars
set(gca,'FontSize',11,'FontName','Arial');
xfac = 1;
yfac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac yfac]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac yfac]);
patch([condvec fliplr(condvec)],[VSO1_mean+VSO1_ste fliplr(VSO1_mean-VSO1_ste)],[0.7 0.7 0.7], 'Linestyle','None'); hold on;
errorbar(condvec, datacond_mean, datacond_ste, 'ok');
xlabel('Condition'); ylabel('Proportion correct');axis([0.8 4.2 0.4 1.0])
set(gca,'YTick',0.4:.2:1.0)
set(gca,'XTick', 1:1:4)

%VEO
figure; %graphs monkey performance across Conditions with error bars
set(gca,'FontSize',11,'FontName','Arial');
xfac = 1;
yfac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac yfac]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac yfac]);
patch([condvec fliplr(condvec)],[VEO1_mean+VEO1_ste fliplr(VEO1_mean-VEO1_ste)],[0.7 0.7 0.7], 'Linestyle','None'); hold on;
errorbar(condvec, datacond_mean, datacond_ste, 'ok');
xlabel('Condition'); ylabel('Proportion correct');axis([0.8 4.2 0.4 1.0])
set(gca,'YTick',0.4:.2:1.0)
set(gca,'XTick', 1:1:4)

%VVO
figure; %graphs monkey performance across Conditions with error bars
set(gca,'FontSize',11,'FontName','Arial');
xfac = 1;
yfac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac yfac]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac yfac]);
patch([condvec fliplr(condvec)],[VVO1_mean+VVO1_ste fliplr(VVO1_mean-VVO1_ste)],[0.7 0.7 0.7], 'Linestyle','None'); hold on;
errorbar(condvec, datacond_mean, datacond_ste, 'ok');
xlabel('Condition'); ylabel('Proportion correct');axis([0.8 4.2 0.4 1.0])
set(gca,'YTick',0.4:.2:1.0)
set(gca,'XTick', 1:1:4)

%% Change magnitude plot
co = [0 0 1;
      0 0.5 0;
      1 0 0;
      0 0.75 0.75;
      0.75 0 0.75;
      0.75 0.75 0;
      0.25 0.25 0.25];
set(groot,'defaultAxesColorOrder',co)
colorvec = get(gca, 'defaultAxesColorOrder');
colorvec = min(colorvec+.65,1);

%VSO
figure;
set(gca,'FontSize',11,'FontName','Arial');
xfac = 1;
yfac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac yfac]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac yfac]);
for pp = 1:length(condvec)
    patch([deltavec, wrev(deltavec)],[VSO2_mean(pp, :) - VSO2_ste(pp, :), wrev(VSO2_mean(pp, :) + VSO2_ste(pp, :))], colorvec(pp, :),'Linestyle','None');
end
hold on;
set(gca,'YTick',0.2:.2:1.0)
set(gca,'XTick',0:30:90)
errorbar(repmat(deltavec,4,1)', datadelta_mean', datadelta_ste','o'); hold on;
xlabel('Change magnitude'); ylabel('Proportion correct');axis([7 93 0.2 1.0]);
legend(strcat('Condition= ',int2str(condvec')), 4, 'Location', 'Best');

%VEO
figure;
set(gca,'FontSize',11,'FontName','Arial');
xfac = 1;
yfac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac yfac]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac yfac]);
for pp = 1:length(condvec)
    patch([deltavec, wrev(deltavec)],[VEO2_mean(pp, :) - VEO2_ste(pp, :), wrev(VEO2_mean(pp, :) + VEO2_ste(pp, :))], colorvec(pp, :),'Linestyle','None');
end
hold on;
set(gca,'YTick',0.2:.2:1.0)
set(gca,'XTick',0:30:90)
errorbar(repmat(deltavec,4,1)', datadelta_mean', datadelta_ste','o'); hold on;
xlabel('Change magnitude'); ylabel('Proportion correct');axis([7 93 0.2 1.0]);
legend(strcat('Condition= ',int2str(condvec')), 4, 'Location', 'Best');

%VVO
figure;
set(gca,'FontSize',11,'FontName','Arial');
xfac = 1;
yfac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac yfac]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac yfac]);
for pp = 1:length(condvec)
    patch([deltavec, wrev(deltavec)],[VVO2_mean(pp, :) - VVO2_ste(pp, :), wrev(VVO2_mean(pp, :) + VVO2_ste(pp, :))], colorvec(pp, :),'Linestyle','None');
end
hold on;
set(gca,'YTick',0.2:.2:1.0)
set(gca,'XTick',0:30:90)
errorbar(repmat(deltavec,4,1)', datadelta_mean', datadelta_ste','o'); hold on;
xlabel('Change magnitude'); ylabel('Proportion correct');axis([7 93 0.2 1.0]);
legend(strcat('Condition= ',int2str(condvec')), 4, 'Location', 'Best');


% BMC plot
figure
set(gca,'FontSize',11,'FontName','Arial');
xfac = 1;
yfac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac yfac]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac yfac]);

mnames = {'VSO','VVO'};
bar(1:2,BMC_meandiff,'k');
hold on
errorbar(1:2,BMC_meandiff,BMC_diffste,'k','LineStyle','none');
set(gca,'XTick',1:2,'XTickLabel',mnames);
set(gca,'YTick',-200:25:0);
ylabel('Model log likelihood relative to VEO');


