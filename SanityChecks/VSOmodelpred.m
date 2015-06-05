clear all; close all;

tic

deltavec = 10:10:90; % change magnitude in degrees
condvec = 1:4; %type of trial condition
Ntrials = 10000; % number of trials per condition
nsteps = 3;

%Parameter ranges for monkey data
% Jbarhvec = linspace(0,100,nsteps);
% Jbarlvec = linspace(0,100,nsteps);
% tauvec = linspace(0.1,100,nsteps);

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
                    
                    %same encoding as other two models
                    %condition 1: both stimuli are HR
                    if condind == 1
                        Jvech = gamrnd(Jbarh/tau, tau, Ntrials, 2);
                        % J to k transformation
                        kappavec = interp1(J_vec,k_vec,Jvech);
                        x = circ_vmrnd(zeros(Ntrials,2),kappavec);
                        
                    %condition 2: changed stimulus is HR and unchanged stimulus is LR
                    elseif condind == 2
                        Jvech = gamrnd(Jbarh/tau, tau, Ntrials, 1);
                        kappavec(:,1) = interp1(J_vec,k_vec,Jvech);
                        x(:,1) = circ_vmrnd(zeros(Ntrials,1),kappavec(:,1));
                        Jvecl = gamrnd(Jbarl/tau, tau, Ntrials, 1);
                        kappavec(:,2) = interp1(J_vec,k_vec,Jvecl);
                        x(:,2) = circ_vmrnd(zeros(Ntrials,1),kappavec(:,2));
                        
                    %condition 3: changed stimulus is LR and unchanged stimulus is HR
                    elseif condind == 3
                        Jvecl = gamrnd(Jbarl/tau, tau, Ntrials, 1);
                        kappavec(:,1) = interp1(J_vec,k_vec,Jvecl);
                        x(:,1) = circ_vmrnd(zeros(Ntrials,1),kappavec(:,1));
                        Jvech = gamrnd(Jbarh/tau, tau, Ntrials, 1);
                        kappavec(:,2) = interp1(J_vec,k_vec,Jvech);
                        x(:,2) = circ_vmrnd(zeros(Ntrials,1),kappavec(:,2));
                        
                    %condition 4: both stimuli are LR
                    elseif condind == 4
                        Jvecl = gamrnd(Jbarl/tau, tau, Ntrials, 2);
                        kappavec = interp1(J_vec,k_vec,Jvecl);
                        x = circ_vmrnd(zeros(Ntrials,2),kappavec);
                    end
                    
                    % Decision rule
                    d = cos (x(:,2)) - (cos(x(:,1)-delta)); %use a single value of kappa so they cancel out
                    
                    perfmodel(condind,deltaind,Jbarhind,Jbarlind,tauind) = (sum(d>0) + 0.5 * sum(d==0))/Ntrials;
                    
                end
            end
        end
        
    end
end

toc
perfmodel(perfmodel==1) = 1-1/Ntrials;
perfmodel(perfmodel==0) = 1/Ntrials;
perfmodel_VSO = perfmodel;

save VSOmodelpred_H %perfmodel_VSO Jbarhvec Jbarlvec tauvec nsteps