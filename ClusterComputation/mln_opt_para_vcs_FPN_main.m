function mln_opt_para_vcs_FPN_main(dirname,filedata,wins)
%%@huifang Aug 23, varying cs, specify the best AUC

dirname='/Users/huifangwang/MULANa/humanCon/HC84S100p32/';
filedata = 'erpN84S1N5100.mat';
wins='1500';


wins=str2double(wins);
winssub=[dirname,'Win', num2str(wins)];
if exist(winssub,'dir')
    
traindatafile = [dirname,'Win', num2str(wins), '/ToutResults/Tout_', num2str(wins),'_',filedata];
else
 traindatafile = [dirname, '/ToutResults/Tout_', num2str(wins),'_',filedata];
   
end
para.name = {'th11','th12','th21','th22','r1','r2','r3','r4','r5','thm'};
para.p0=[0.5,0.5,0.6,0.6,0,0,0,0,1,0.5];
para.rangeMin=[0.2,0.2,0.3,0.3,0,0,0,0,0.8,0.5];
para.rangeMax=[0.9,0.9,0.9,0.9,0.5,0.5,0.5,0.5,1,0.9];
fitness={'FPN'};

%methodlog={'BCorrU','COH1','BCorrD','PDC'};
methodgroup={{'BTEU','BCorrU','hmvar','ffDTF'},{'BCorrU','BTEU','BTED','ffDTF'},...,
    {'BCorrU','PCorrU','hmvar','ffDTF'},{'BCorrU','PCorrU','BTED','hmvar'},{'BCorrU','COH1','BCorrD','PDC'}};

for ibg=[2,5]%1:length(methodgroup)
para.fitness=fitness{1};
methodlog=methodgroup{ibg};
savefilefit = [dirname,'Win', num2str(wins), '/ToutResults/Fit_',para.fitness, 'G',num2str(ibg),'_',num2str(wins),'_',filedata];
mln_opt_para_vcs_fpn(traindatafile,para,methodlog,savefilefit);
end

