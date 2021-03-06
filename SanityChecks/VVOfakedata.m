function fakedata = VVOfakedata(Jbarh, Jbarl, tau)

% This code creates a fake data set for the VPF model

Ntrials = 2500; % number of trials per condition
condvec = 1:4;
deltavec = 10:10:90; % change magnitude in degrees

% interpolation of kappa over J
upper_bound = 3000;
k_vec = linspace(0,upper_bound,1e4);
J_vec = k_vec .* besseli(1,k_vec,1)./besseli(0,k_vec,1);

fakedata = [];
for condind = 1:length(condvec)
    cond = condvec(condind);
    
    fakedataN = zeros(Ntrials,3);
    % Generative model: simulating observer's measurements
    
    fakedataN(:,1) = cond;
    delta = deltavec(randi(length(deltavec),[1 Ntrials]))'* 2*pi/180; % randomly generating change magnitude from deltavec
    %%
    
    %condition 1: both stimuli are HR
    if condind == 1
        Jvech = gamrnd(Jbarh/tau, tau, Ntrials, 2);
        kappavec = interp1(J_vec,k_vec,Jvech);%draw both kappavecs from Jbar high for encoding and decision rule
        x = circ_vmrnd(zeros(Ntrials,2),kappavec);
        
        %condition 2: changed stimulus is HR and unchanged stimulus is LR
    elseif condind == 2
        Jvech = gamrnd(Jbarh/tau, tau, Ntrials, 1);
        kappavec(:,1) = interp1(J_vec,k_vec,Jvech);% draw HR kappavec from jbar high for encoding and decision rule
        x(:,1) = circ_vmrnd(zeros(Ntrials,1),kappavec(:,1));
        Jvecl = gamrnd(Jbarl/tau, tau, Ntrials, 1);
        kappavec(:,2) = interp1(J_vec,k_vec,Jvecl);% draw LR kappavec from jbar low for encoding and decision rule
        x(:,2) = circ_vmrnd(zeros(Ntrials,1),kappavec(:,2));
        
        %condition 3: changed stimulus is LR and unchanged stimulus is HR
    elseif condind == 3
        Jvecl = gamrnd(Jbarl/tau, tau, Ntrials, 1);
        kappavec(:,1) = interp1(J_vec,k_vec,Jvecl);% draw LR kappavec from jbar low for encoding and decision rule
        x(:,1) = circ_vmrnd(zeros(Ntrials,1),kappavec(:,1));
        Jvech = gamrnd(Jbarh/tau, tau, Ntrials, 1);
        kappavec(:,2) = interp1(J_vec,k_vec,Jvech);% draw HR kappavec from jbar high for encoding and decision rule
        x(:,2) = circ_vmrnd(zeros(Ntrials,1),kappavec(:,2));
        
        %condition 4: both stimuli are LR
    elseif condind == 4
        Jvecl = gamrnd(Jbarl/tau, tau, Ntrials, 2);
        kappavec = interp1(J_vec,k_vec,Jvecl);% draw both kappavecs from jbar low for encoding and decision rule
        x = circ_vmrnd(zeros(Ntrials,2),kappavec);
    end
    
    % Decision rule
    d = -log (besseli(0, kappavec(:,2)))+ kappavec(:,2) .* cos (x(:,2)) + log (besseli(0, kappavec(:,1)))- kappavec(:,1) .* cos(x(:,1)-delta);
    fakedataN(:,2) = round(delta/(2*pi/180));
    fakedataN(:,3) = (d>0) + (d==0).* round(rand(Ntrials,1));
    
    
    fakedata = [fakedata ; fakedataN];
end

%% calculating and plotting fakedata performance across set sizes and change magnitudes

perfdata1 = zeros(1,length(condvec));
perfdata2  = zeros(length(condvec),length(deltavec));

conddata =fakedata(:,1);
deltadata = fakedata(:,2);
corrdata = fakedata(:,3);

for condind = 1:length(condvec)
    cond = condvec(condind);
    perfdata1(condind) = mean(corrdata(conddata == cond)); % perfdata1 calculates performance for different set sizes
    
    for deltaind = 1:length(deltavec)
        delta = deltavec(deltaind);
        perfdata2(condind,deltaind) = mean(corrdata(conddata == cond & deltadata == delta)); % perfdata2 calculates performance for different set sizes and change in magnitude
    end
end


datacond_mean = perfdata1;
datadelta_mean = perfdata2;

figure;
plot(condvec, datacond_mean);
xlabel('Set size'); ylabel('Proportion correct');axis([0.8 4.2 0.5 1.0])
set(gca,'YTick',0.4:.1:1.0)
set(gca,'XTick', 1:1:4)

colorvec = get(gca, 'ColorOrder');
colorvec = min(colorvec+.65,1);

figure;
set(gca,'YTick',0.4:.1:1)
set(gca,'XTick',10:10:90)
plot(repmat(deltavec,4,1)', datadelta_mean','-o'); hold on;
xlabel('Change magnitude'); ylabel('Proportion correct');axis([10 90 0.4 1]);
legend(strcat('cond= ',int2str(condvec')), 4);

save fakedataVVO.mat perfdata1 perfdata2
