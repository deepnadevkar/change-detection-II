function VEOmodelpred(Jbarhind)

tic
deltavec = 10:10:90; % change magnitude in degrees
condvec = 1:4; %type of trial condition
Ntrials = 10000; % number of trials per condition
nsteps = 100;

%Parameter ranges
Jbarhvec = linspace(0,100,nsteps);
Jbarlvec = linspace(0,100,nsteps);
tauvec = linspace(0.1,100,nsteps);

% interpolation of kappa over J
upper_bound = 3000;
k_vec = linspace(0,upper_bound,1e4);
J_vec = k_vec .* besseli(1,k_vec,1)./besseli(0,k_vec,1);

%% Calculate model predictions and fitting model

perfmodel = NaN(length(condvec),length(deltavec),length(Jbarlvec),length(tauvec));

Jbarh = Jbarhvec(Jbarhind);
kappah = interp1(J_vec, k_vec, Jbarh);

for Jbarlind = 1:length(Jbarlvec)
    Jbarl = Jbarlvec(Jbarlind);
    
    kappal = interp1(J_vec, k_vec, Jbarl);
    
    for tauind = 1:length(tauvec)
        tau = tauvec(tauind);
        
        for deltaind = 1: length(deltavec)
            delta = deltavec(deltaind) * 2* pi/180;
            
            for condind = 1:length(condvec)
                
                if condind == 1 %condition 1: both stimuli are HR
                    
                    Jvech = gamrnd(Jbarh/tau, tau, Ntrials, 2);
                    kappavec = interp1(J_vec,k_vec,Jvech); %draw both kappavecs from Jbar high for encoding
                    x = circ_vmrnd(zeros(Ntrials,2),kappavec);
                    
                    kappa1 = kappah; %draw both kappas from Jbar high for decision rule
                    kappa2 = kappa1;
                    
                elseif condind == 2  %condition 2: changed stimulus is HR and unchanged stimulus is LR
                    Jvech = gamrnd(Jbarh/tau, tau, Ntrials, 1);
                    kappavec(:,1) = interp1(J_vec,k_vec,Jvech);  %draw kappavec for HR from Jbar high for encoding
                    x(:,1) = circ_vmrnd(zeros(Ntrials,1),kappavec(:,1));
                    Jvecl = gamrnd(Jbarl/tau, tau, Ntrials, 1);
                    kappavec(:,2) = interp1(J_vec,k_vec,Jvecl); %draw kappavec for LR from Jbar low for encoding
                    x(:,2) = circ_vmrnd(zeros(Ntrials,1),kappavec(:,2));
                    
                    kappa1 = kappah; %draw kappa high for decision rule
                    kappa2 = kappal; %draw kappa low for decision rule                    
                    
                elseif condind == 3 %condition 3: changed stimulus is LR and unchanged stimulus is HR
                    Jvecl = gamrnd(Jbarl/tau, tau, Ntrials, 1);
                    kappavec(:,1) = interp1(J_vec,k_vec,Jvecl); %draw kappavec for LR from Jbar low for encoding
                    x(:,1) = circ_vmrnd(zeros(Ntrials,1),kappavec(:,1));
                    Jvech = gamrnd(Jbarh/tau, tau, Ntrials, 1);
                    kappavec(:,2) = interp1(J_vec,k_vec,Jvech);%draw kappavec for HR from Jbar high for encoding
                    x(:,2) = circ_vmrnd(zeros(Ntrials,1),kappavec(:,2));
                    
                    kappa1 = kappal; %draw kappa low for decision rule
                    kappa2 = kappah; %draw kappa high for decision rule
                    
                    
                elseif condind == 4 %condition 4: both stimuli are LR
                    Jvecl = gamrnd(Jbarl/tau, tau, Ntrials, 2);
                    kappavec = interp1(J_vec,k_vec,Jvecl);
                    x = circ_vmrnd(zeros(Ntrials,2),kappavec);
                    
                    kappa1 = kappal; %draw both kappas from Jbar low for decision rule
                    kappa2 = kappa1;
                end
                
                % Decision rule
                d = -log (besseli(0, kappa2))+ kappa2 .* cos (x(:,2)) + log (besseli(0, kappa1))- kappa1 .* cos(x(:,1)-delta);
                
                perfmodel(condind,deltaind,Jbarlind,tauind) = (sum(d>0) + 0.5 * sum(d==0))/Ntrials;
                
            end
        end
    end
end

toc
perfmodel(perfmodel==1) = 1-1/Ntrials;
perfmodel(perfmodel==0) = 1/Ntrials;
perfmodel_VEO = perfmodel;

filename = strcat('VEOmodelpred',num2str(Jbarhind))
save(filename, 'perfmodel_VEO','Jbarh','Jbarlvec','tauvec')