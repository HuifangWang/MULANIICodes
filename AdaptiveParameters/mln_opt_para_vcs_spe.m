function mln_opt_para_vcs_spe(traindatafile,para,methodlog,savefilefit,bestvcs)
% Huifang Wang, May/28/2015;  find traindata from traindatafile
MatA = load(traindatafile,methodlog{:},'Connectivity');

InputMM=mln_GenerateTDM_win(MatA,methodlog,100,40);

iInput=InputMM.values(:,:,1);
fMat = MatA.Connectivity;
kforMF=5;
PopulationSize=16;
options.PopulationSize = PopulationSize;
GenomeLength = length(para.p0);
options.InitialPopulation=para.p0;
options.PopulationType = 'doubleVector';
options.PopInitRangeMin = para.rangeMin;
options.PopInitRangeMax = para.rangeMax;
options.EliteCount=2;
options.CrossoverFraction=0.8;
options.MutationGens = 2;
para.PopulationSize=options.PopulationSize;
fitness=para.fitness;
minfMat=min(fMat(fMat>=0.1));
maxfMat=max(max(fMat));
vcs=minfMat:0.1:maxfMat;
indvcsbest=mln_npindex(vcs,bestvcs);
state.vcs=vcs;
% initialize the first population
Population = mln_gacreationuniform(GenomeLength,options);
para.pc=Population;
state.Generation=1;
exitFlag=0;
while ~exitFlag
Population=para.pc;
[igAUCvcs,igthetaAUC,Population,Matall] = mln_GA_forward(iInput,Population,para,kforMF,methodlog,fMat,vcs);
%use scorepfit or scoreAUC
switch fitness
    case 'AUCOS'
    scorepfit=max(igAUCvcs,[],1);
    case 'AUCspecs'
        igscorepfit=igAUCvcs(indvcsbest,:);
        scorepfit=mean(igscorepfit,1);
        
end
 [state.bestscorepfit(state.Generation),indbest] =max(scorepfit);
 state.bestMat(:,:,state.Generation)=Matall(:,:,indbest);
 state.bestCh(state.Generation,:) = Population(indbest,:);
state.AUCvcs(state.Generation,:,:) = igAUCvcs;
state.thetaAUC(state.Generation,:,:) = igthetaAUC;

state.pc= Population;
%state=mln_plot_GA(state);

exitFlag = mln_GA_isTimetoStop(state);
    if ~exitFlag 
        
% update GA population
    [para.pc,state] = mln_stepGA(scorepfit,Population,options,state,GenomeLength);
    state.Generation=state.Generation+1;
    end
    
end

save(savefilefit,'state','options');






function exitFlag = mln_GA_isTimetoStop(state)
exitFlag=0;
if state.bestscorepfit(state.Generation)>0.999
    exitFlag=1;
    return;
end
if state.Generation >3
if var(state.bestscorepfit(state.Generation-3:state.Generation))<1e-9
    exitFlag=1;
    return;
end
end

if state.Generation > 40
exitFlag=1;
    return;
end
function Input=mln_GenerateTDM_win(MatA,methods_used,nboot,nboot2)

Nchan=size(MatA.(methods_used{1}),1);

Nlink=Nchan*(Nchan-1);
Nmethod=length(methods_used);
Input.values=NaN(Nlink,Nmethod,nboot2);
Input.methods=methods_used;

nbins=100;

for imethod=1:Nmethod
    Mat1=MatA.(methods_used{imethod});
    [Nnode,~,Nwin]=size(Mat1);
    isizeboot=Nlink;
    
    MatBoot=NaN( isizeboot*nboot,1);
    
    iMat=abs(mean(Mat1,3));
    iMat(1:size(iMat,1)+1:end) = NaN;
    a=reshape(iMat,Nnode*Nnode,1);
    MatBoot(1:isizeboot,1)=a(~isnan(a));
    
    for iboot=2:nboot
        iBoot_ind=randi(Nwin,[1,Nwin]);
        iBoot_Mat=Mat1(:,:,iBoot_ind);
        iMat=abs(mean(iBoot_Mat,3));
        iMat(1:size(iMat,1)+1:end) = NaN;
        a=reshape(iMat,Nnode*Nnode,1);
        MatBoot((iboot-1)*isizeboot+1:iboot*isizeboot,1)=a(~isnan(a));
    end
    xi=linspace(min(MatBoot),max(MatBoot),nbins);
    n_elements = histc(MatBoot,xi);
    c_elements = cumsum(n_elements);
    c_elements=c_elements./max(c_elements);
    %% seperate MatBoot to MBLink=M(Nlink,)
    %% the means for all windows
    MBLink=reshape(MatBoot,Nlink,nboot);
    for iboot=1:nboot2
        for ipairs=1:Nlink
            Input.values(ipairs,imethod,iboot)=c_elements(find(xi>=MBLink(ipairs,iboot),1));
        end
    end

end

function [AUCvcs,thetaAUC,Population,Matall] = mln_GA_forward(iInput,Population,para,kforMF,methodlog,fMat,vcs)

% indthf=strcmp(para.name,'thf');
% rangethf = [para.rangeMin(indthf),para.rangeMax(indthf)];
%calculate within each Population
%pfit=zeros(1,para.PopulationSize);
%AUCfit=zeros(1,para.PopulationSize);
Matall=zeros([size(fMat),size(Population,1)]);


AUCvcs=NaN(length(vcs),para.PopulationSize);
thetaAUC=NaN(length(vcs),para.PopulationSize);

for iPop = 1: para.PopulationSize
para.pc = Population(iPop,:);    
[iMat,~] = mln_algorithm_para(iInput,para,kforMF,methodlog);

for ivcs=1:length(vcs)
[AUCvcs(ivcs,iPop),thetaAUC(ivcs,iPop)]=mln_AUC(iMat,fMat.*(fMat>(vcs(ivcs)-0.001)),0,1);
end

% iMatRef=iMat.*(iMat>=para.pc(indthf));   
% ipfit=mln_Sim_Cos(iMatRef,fMat);
% if iAUCfit>0.9 && ipfit<0.8
%     [thf,ipfit]=mln_updata_thf(iMat,fMat,rangethf,iAUCfit);
%     Population(iPop,indthf)=thf;
% end
%pfit(iPop)=ipfit;
%AUCfit(iPop)=iAUCfit;
Matall(:,:,iPop)=iMat;
end


function state=mln_plot_GA(state)

if state.Generation==1
   state.figplot=figure; 
end

np=size(state.AUCvcs,3);

set(gca,'Parent',state.figplot)
ax=gca;
hold on
for ipop=1:np
plot(state.vcs,state.AUCvcs(state.Generation,:,ipop),'bo','MarkerFaceColor','b');
plot(state.vcs,state.AUCvcs(state.Generation,:,ipop),'b-');

plot(state.vcs+0.2,state.thetaAUC(state.Generation,:,ipop),'m*','MarkerFaceColor','b');
plot(state.vcs+0.2,state.thetaAUC(state.Generation,:,ipop),'m-');
end
