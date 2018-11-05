function mln_opt_para_vcs_spe_main(dirname,filedata,wins)
%%@huifang Aug 23, varying cs, specify the best AUC

dirname='/Users/huifangwang/MULANa/mlnData/N20CS1/';
filedata = 'nmmN20L10CS10100S7N5100.mat';
wins=1000;
%wins=str2double(wins);
traindatafile = [dirname,'Win', num2str(wins), '/ToutResults/Tout_', num2str(wins),'_',filedata];
para.name = {'th11','th12','th21','th22','r1','r2','r3','r4','r5','thm'};
para.p0=[0.5,0.5,0.6,0.6,0,0,0,0,1,0.5];
para.rangeMin=[0.2,0.2,0.3,0.3,0,0,0,0,0.8,0.5];
para.rangeMax=[0.9,0.9,0.9,0.9,0.5,0.5,0.5,0.5,1,0.9];
fitness={'AUCspecs'};
bestvcs=[0.5,0.6,0.7];

methodlog={'BCorrU','COH1','BCorrD','PDC'};

for ifit = 1: length(fitness)
    para.fitness=fitness{ifit};
savefilefit = [dirname,'Win', num2str(wins), '/ToutResults/Fit_',para.fitness, num2str(wins),'_',filedata];
mln_opt_para_vcs_spe(traindatafile,para,methodlog,savefilefit,bestvcs);
end
