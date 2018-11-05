function mln_opt_para_main(dirname,filedata,wins)
%dirname='/Users/huifangwang/MULANa/mlnData/N10CS8/';
%filedata = 'nmmN10L5CS80100S7N5100.mat';
%wins=1000;
%wins=str2double(wins);
traindatafile = [dirname,'Win', num2str(wins), '/ToutResults/Tout_', num2str(wins),'_',filedata];
para.name = {'th11','th12','th21','th22','r1','r2','r3','r4','r5','thm','thf'};
para.p0=[0.5,0.5,0.6,0.6,0,0,0,0,1,0.5,0.999];
para.rangeMin=[0.2,0.2,0.3,0.3,0,0,0,0,0.8,0.5,0.8];
para.rangeMax=[0.9,0.9,0.9,0.9,0.5,0.5,0.5,0.5,1,0.9,0.999];
fitness={'AUC','COS'};


methodlog={'BCorrU','COH1','BCorrD','PDC'};

for ifit = 1: length(fitness)
    para.fitness=fitness{ifit};
savefilefit = [dirname,'Win', num2str(wins), '/ToutResults/Fit_',para.fitness, num2str(wins),'_',filedata];
mln_opt_para(traindatafile,para,methodlog,savefilefit);
end
