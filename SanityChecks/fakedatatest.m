function [data1, data2, VSO1, VSO2, VEO1, VEO2, VVO1, VVO2, pars, parsest_VSO, parsest_VEO, parsest_VVO, LL, BIC, AIC, AICc, BMC] = fakedatatest(generatingmodel)

%% GENERATING DATA

condvec = 1:4;
deltavec = 10:10:90; % change magnitude in degrees
nsteps = 100;

switch generatingmodel
    case 1
        Jbarhvec = linspace(0,100,nsteps);
        Jbarlvec = linspace(0,100,nsteps);
        tauvec = linspace(0,100,nsteps);
        Jbarh = Jbarhvec(randi(length(Jbarhvec),1));   
        Jbarl = Jbarlvec(randi(length(Jbarlvec),1)); 
        tau = tauvec(randi(length(tauvec),1));
        pars = [Jbarh Jbarl tau];
        fakedata = VSOfakedata(Jbarh,Jbarl,tau);
    case 2
        Jbarhvec = linspace(0,100,nsteps);
        Jbarlvec = linspace(0,100,nsteps);
        tauvec = linspace(0,100,nsteps);
        Jbarh = Jbarhvec(randi(length(Jbarhvec),1));   
        Jbarl = Jbarlvec(randi(length(Jbarlvec),1)); 
        tau = tauvec(randi(length(tauvec),1));
        pars = [Jbarh Jbarl tau];
        fakedata = VEOfakedata(Jbarh,Jbarl,tau);
    case 3
        Jbarhvec = linspace(0,100,nsteps);
        Jbarlvec = linspace(0,100,nsteps);
        tauvec = linspace(0,100,nsteps);
        Jbarh = Jbarhvec(randi(length(Jbarhvec),1));   
        Jbarl = Jbarlvec(randi(length(Jbarlvec),1)); 
        tau = tauvec(randi(length(tauvec),1));
        pars = [Jbarh Jbarl tau];
        fakedata = VVOfakedata(Jbarh,Jbarl,tau);
end
   
%% calculating fakedata performance across conditions and change magnitudes
conddata = fakedata(:,1);
deltadata = fakedata(:,2);
corrdata = fakedata(:,3);

data1 = zeros(length(condvec),1);
data2 = zeros(length(condvec),length(deltavec));
data1_ste = zeros(length(condvec),1);
data2_ste = zeros(length(condvec),length(deltavec));

for Nind = 1:length(condvec)
    N = condvec(Nind);
    data1(Nind) = mean(corrdata(conddata == N)); % data1 calculates ormance for different set sizes
    data1_ste (Nind) = std(corrdata(conddata == N));
    for deltaind = 1:length(deltavec)
        delta = deltavec(deltaind);
        data2(Nind,deltaind) = mean(corrdata(conddata == N & deltadata == delta)); % data2 calculates ormance for different set sizes and change in magnitude
        data2_ste (Nind, deltaind) = std(corrdata(conddata == N & deltadata == delta));
    end
end

%% TESTING MODELS ON THE GENERATED DATA

% parsest_VSO = [Jbarhest Jbarlest tauest];
% parsest_VEO = [Jbarhest Jbarlest tauest];
% parsest_VEO = [Jbarhest Jbarlest tauest];

[VSO1, VSO2, VEO1, VEO2, VVO1, VVO2, parsest_VSO, parsest_VEO, parsest_VVO, LL, BMC, BIC, AIC, AICc, Rs_VSO, Rs_VEO, Rs_VVO] = modelfitting(fakedata);


