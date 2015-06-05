clear all; close all;

tic

deltavec = 10:10:90; % change magnitude in degrees
condvec = 1:4; %type of trial condition
Ntrials = 10000; % number of trials per condition
nsteps = 3;

%Parameter ranges for monkey data
% Jbarhvec = linspace(0,30,nsteps);
% Jbarlvec = linspace(0,30,nsteps);
% tauvec = linspace(0.1,30,nsteps);

%Parameter ranges for human data
Jbarhvec = linspace(0,100,nsteps);
Jbarlvec = linspace(0,100,nsteps);
tauvec = linspace(0.1,100,nsteps);

% interpolation of kappa over J
upper_bound = 3000;
k_vec = linspace(0,upper_bound,1e4);
J_vec = k_vec .* besseli(1,k_vec,1)./besseli(0,k_vec,1);

%% Calculate model predictions and fitting model

perfmodel = zeros(length(condvec),length(deltavec),length(Jbarhvec),length(Jbarlvec),length(tauvec));

for Jbarhind = 1:length(Jbarhvec)
    Jbarh = Jbarhvec(Jbarhind);
    
    for Jbarlind = 1:length(Jbarlvec)
        Jbarl = Jbarlvec(Jbarlind);
        
        for tauind = 1:length(tauvec)
            tau = tauvec(tauind);
            
            for deltaind = 1: length(deltavec)
                delta = deltavec(deltaind) * 2* pi/180;
                
                for condind = 1:length(condvec)
                    cond = condvec(condind);
                                    
                    %condition 1: both stimuli are HR
                    if condind == 1
                        Jvech = gamrnd(Jbarh/tau, tau, Ntrials, 2);
                        kappavec = interp1(J_vec,k_vec,Jvech); %draw both kappavecs from Jbar high for encoding
                        x = circ_vmrnd(zeros(Ntrials,2),kappavec);
                        
                        kappa1 = interp1(J_vec, k_vec, Jbarh); %draw both kappas from Jbar high for decision rule
                        kappa2 = kappa1;
                        
                    %condition 2: changed stimulus is HR and unchanged stimulus is LR
                    elseif condind == 2
                        Jvech = gamrnd(Jbarh/tau, tau, Ntrials, 1);
                        kappavec(:,1) = interp1(J_vec,k_vec,Jvech);  %draw kappavec for HR from Jbar high for encoding
                        x(:,1) = circ_vmrnd(zeros(Ntrials,1),kappavec(:,1));
                        Jvecl = gamrnd(Jbarl/tau, tau, Ntrials, 1);
                        kappavec(:,2) = interp1(J_vec,k_vec,Jvecl); %draw kappavec for LR from Jbar low for encoding
                        x(:,2) = circ_vmrnd(zeros(Ntrials,1),kappavec(:,2));
                        
                        kappa1 = interp1(J_vec, k_vec, Jbarh); %draw kappa high for decision rule
                        kappa2 = interp1(J_vec, k_vec, Jbarl); %draw kappa low for decision rule

                    %condition 3: changed stimulus is LR and unchanged stimulus is HR
                    elseif condind == 3
                        Jvecl = gamrnd(Jbarl/tau, tau, Ntrials, 1);
                        kappavec(:,1) = interp1(J_vec,k_vec,Jvecl); %draw kappavec for LR from Jbar low for encoding
                        x(:,1) = circ_vmrnd(zeros(Ntrials,1),kappavec(:,1));
                        Jvech = gamrnd(Jbarh/tau, tau, Ntrials, 1);
                        kappavec(:,2) = interp1(J_vec,k_vec,Jvech);%draw kappavec for HR from Jbar high for encoding
                        x(:,2) = circ_vmrnd(zeros(Ntrials,1),kappavec(:,2));
                                                
                        kappa1 = interp1(J_vec, k_vec, Jbarl);%draw kappa low for decision rule
                        kappa2 = interp1(J_vec, k_vec, Jbarh);%draw kappa high for decision rule

                    %condition 4: both stimuli are LR
                    elseif condind == 4
                        Jvecl = gamrnd(Jbarl/tau, tau, Ntrials, 2);
                        kappavec = interp1(J_vec,k_vec,Jvecl);
                        x = circ_vmrnd(zeros(Ntrials,2),kappavec);
                                                
                        kappa1 = interp1(J_vec, k_vec, Jbarl);%draw both kappas from Jbar low for decision rule
                        kappa2 = kappa1;
                    end
                    
                    % Decision rule
                    d = -log (besseli(0, kappa2))+ kappa2 .* cos (x(:,2)) + log (besseli(0, kappa1))- kappa1 .* cos(x(:,1)-delta);
                    
                    perfmodel(condind,deltaind,Jbarhind,Jbarlind,tauind) = (sum(d>0) + 0.5 * sum(d==0))/Ntrials;
                    
                end
            end
        end
            end
end

toc
perfmodel(perfmodel==1) = 1-1/Ntrials;
perfmodel(perfmodel==0) = 1/Ntrials;
perfmodel_VEO = perfmodel;

save VEOmodelpred_H %perfmodel_VEO Jbarhvec Jbarlvec tauvec nsteps