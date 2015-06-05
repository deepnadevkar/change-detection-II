clear all; close all;
condvec = 1:4;
deltavec = 10:10:90;

generatingmodel = 1;
tic
[data1, data2, VSO1, VSO2, VEO1, VEO2, VVO1, VVO2, pars, parsest_VSO, parsest_VEO, parsest_VVO, LL, BIC, AIC, AICc, BMC] = fakedatatest(generatingmodel);
% Parameter recovery & difference in criterions from winning model 
switch generatingmodel
    case 1
        parsest = parsest_VSO;
        
    case 2
        parsest = parsest_VEO;
        
    case 3
        parsest = parsest_VVO;
          
end

[pars; parsest]

% VSO model fits
figure; plot(condvec,data1,'-o')
hold on; plot(condvec,VSO1, '--');
xlabel('Set size'); ylabel('Proportion correct');axis([0.8 4.2 0.3 1.1])
set(gca,'YTick',0.3:.1:1.0)
set(gca,'XTick', 1:1:4)
title_VSO_ = title('VSO model fit');
set(title_VSO_,'String','VSO model fit')

figure;
plot(deltavec,data2'); hold on;
plot(deltavec, VSO2','--');
xlabel('Change magnitude'); ylabel('Proportion correct');
ylim([0.3 1.3])
legend(strcat('Condition =',num2str(condvec')),'Location','Best')
title_VSO_ = title('VSO model fit');
set(title_VSO_,'String','VSOmodel fit')

% VEO model fits
figure; plot(condvec,data1,'-o')
hold on; plot(condvec,VEO1, '--');
xlabel('Set size'); ylabel('Proportion correct');axis([0.8 4.2 0.3 1.1])
set(gca,'YTick',0.3:.1:1.0)
set(gca,'XTick', 1:1:4)
title_VEO_ = title('VEO model fit');
set(title_VEO_,'String','VEO model fit')

figure;
plot(deltavec,data2'); hold on;
plot(deltavec, VEO2','--');
xlabel('Change magnitude'); ylabel('Proportion correct');
ylim([0.3 1.3])
legend(strcat('Condition =',num2str(condvec')),'Location','Best')
title_VEO_ = title('VEO model fit');
set(title_VEO_,'String','VEOmodel fit')

% VVO model fits
figure; plot(condvec,data1,'-o')
hold on; plot(condvec,VVO1, '--');
xlabel('Set size'); ylabel('Proportion correct');axis([0.8 4.2 0.3 1.1])
set(gca,'YTick',0.3:.1:1.0)
set(gca,'XTick', 1:1:4)
title_VVO_ = title('VVO model fit');
set(title_VVO_,'String','VVO model fit')

figure;
plot(deltavec,data2'); hold on;
plot(deltavec, VVO2','--');
xlabel('Change magnitude'); ylabel('Proportion correct');
ylim([0.3 1.3])
legend(strcat('Condition =',num2str(condvec')),'Location','Best')
title_VVO_ = title('VVO model fit');
set(title_VVO_,'String','VVOmodel fit')

toc






