function mln_save_para_gb()

%% save the calculated parameters
% this one is special for only cs is large than 8
vnc=[20,30];
nl=15;
Nstas=50;

methodgroup={{'BTEU','BCorrU','hmvar','ffDTF'},{'BCorrU','BTEU','BTED','ffDTF'},...,
    {'BCorrU','PCorrU','hmvar','ffDTF'},{'BCorrU','PCorrU','BTED','hmvar'},{'BCorrU','COH1','BCorrD','PDC'}};
BG=[2,5];

cs=8;
Nbins=20;
for indnc =1:2
    for indBG =1:length(BG)
    iG=BG(indBG);
    nc=vnc(indnc);
    basicdir=['/Volumes/Lioncub/Users/huifangwang/MULANa/mlnData/N',num2str(nc),'CS',num2str(cs),'S',num2str(Nstas),'/'];
    iwins=1000;
    CSm=100;

    fitness='FPN';
    showBestFPN=[];
    for istru =1:Nstas
        datasets=['nmmN',num2str(nc),'L',num2str(nl),'CS',num2str(cs*10),num2str(CSm),'S',num2str(istru),'N5100.mat'];
        filename=[basicdir,'Win1000','/ToutResults/','Fit_',fitness,'G',num2str(iG),'_',num2str(iwins),'_',datasets];
        load(filename,'state');
        
        if state.bestscorepfit>-0.1
            ind0=find(state.FPN==0);
            showBest=[state.pc(ind0,:),state.theta(end,ind0)'];% please check this 1 here
            if isempty(showBestFPN)
               showBestFPN = showBest;
            else
               showBestFPN = [showBestFPN;showBest];
            end
        end
    end
    
   opt(indBG).BMs=methodgroup{indBG};
   opt(indBG).iG=iG;
   opt(indBG).BestFPNch=showBestFPN;
   opt(indBG).PE=median(showBestFPN,1);
   opt(indBG).var=var(showBestFPN,1)/2;
   opt(indBG).tak=showBestFPN;
   opt(indBG).mean=mean(showBestFPN,1);
   opt(indBG).hp=mln_findpk(showBestFPN,1,Nbins);
   savefilename=[basicdir,'opt_para_N',num2str(nc),'CS',num2str(cs),'L',num2str(nl),'.mat'];
   save(savefilename,'opt')
end
end


function hp=mln_findpk(showBest,ax,Nbins)

[nr,nc]=size(showBest);
if ax==1
    hp=NaN(1,nc);
    for inc=1:nc
        dataset=showBest(:,inc);
        [y,x]=hist(dataset,Nbins);
        [~,indy]=max(y,[],2);
        hp(inc)=x(indy);
    end
    else
    hp=NaN(1,nr);
    for inr=1:nr
        dataset=showBest(inr,:);
        [y,x]=hist(dataset,Nbins);
        [~,indy]=max(y,[],2);
        hp(inr)=x(indy);
    end
end