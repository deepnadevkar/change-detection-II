clear all close all; 
perfmodel = zeros(4,9,100,100,100);

for ii = 1:100
   fname = sprintf('Z:/dthakkar/CDTIICode - fakedatatest_9.1.14/VVO/VVOmodelpredictions/VVOmodelpred%i.mat', ii);
   load(fname);
   perfmodel(:,:,ii,:,:) = perfmodel_VVO;
end

clear perfmodel_VVO
perfmodel_VVO = perfmodel; 

Jbarhvec = linspace(0,100,100);

save modelpred_VVO perfmodel_VVO Jbarhvec Jbarlvec tauvec