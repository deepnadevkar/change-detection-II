clear all close all;

nsteps = 100;
Jbarhvec = linspace(0,100,nsteps);
Jbarlvec = linspace(0,100,nsteps);
tauvec = linspace(0,100,nsteps);

for i = 1:20;
Jbarh(i) = Jbarhvec(randi(length(Jbarhvec),1));
Jbarl(i) = Jbarlvec(randi(length(Jbarlvec),1));
tau(i) = tauvec(randi(length(tauvec),1));

pars_VSO(i,:) = [Jbarh(i) Jbarl(i) tau(i)];
end

for i = 1:20;
Jbarh(i) = Jbarhvec(randi(length(Jbarhvec),1));
Jbarl(i) = Jbarlvec(randi(length(Jbarlvec),1));
tau(i) = tauvec(randi(length(tauvec),1));

pars_VEO(i,:) = [Jbarh(i) Jbarl(i) tau(i)];
end

for i = 1:20;
Jbarh(i) = Jbarhvec(randi(length(Jbarhvec),1));
Jbarl(i) = Jbarlvec(randi(length(Jbarlvec),1));
tau(i) = tauvec(randi(length(tauvec),1));

pars_VVO(i,:) = [Jbarh(i) Jbarl(i) tau(i)];
end
save fakedatapars pars_VSO pars_VEO pars_VVO