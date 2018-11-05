function mln_save_para_gb_vcs_hidden_HH()

%% save the calculated parameters
%vnc=[20,30];
vnc=[30];
ncrit=0;

Nstas=50;
NP=20100;
methodgroup={{'BTEU','BCorrU','hmvar','ffDTF'},{'BCorrU','BTEU','BTED','ffDTF'},...,
    {'BCorrU','PCorrU','hmvar','ffDTF'},{'BCorrU','PCorrU','BTED','hmvar'},{'BCorrU','COH1','BCorrD','PDC'}};
%BG=[1,2,3,4,5];
BG=[2,5];
cs=1;
Nbins=20;

bgv={'12','13','23'};

% 
    
 
for indnc =1:length(vnc)
    for indBG =1:length(BG)
        iG=BG(indBG);
        
        nc=vnc(indnc);
        basicdir=['/Users/huifangwang/MULANa/humanCon/HH',num2str(nc),'S50T/hidden/'];
        iwins=1500;
        
        ds=4;
        fitness='FPN';
        showBestFPN=[];
        for istru =1:50 %istru =[1:4,10:49]%1:Nstas
            datasets=['erpN',num2str(nc),'S',num2str(istru),'N',num2str(NP)];
            filename=[basicdir,'/ToutResults/','Fit_',fitness,'G',num2str(iG),'_',num2str(iwins),'_',datasets];
            if ds>0
                filename=[filename,'ds',num2str(ds)];
            end
            
            for ibgv=1:length(bgv)
                sfilename=[filename,'b',bgv{ibgv}];
                try
                load(sfilename,'state');
                catch
                    continue
                
                end
                if state.bestscorepfit(end)>-0.1-ncrit
                    
                    iFPN=sum(state.FPN(end,:,:),3);
                    ind0=find(iFPN<=ncrit);
                    itheta=squeeze(state.theta(end,:,:));
                    showBest=[state.pc(ind0,:),itheta(ind0,:)];
                    if isempty(showBestFPN)
                        showBestFPN = showBest;
                    else
                        showBestFPN = [showBestFPN;showBest];
                    end
                end
            end
        end
        
        opt(indBG).BMs=methodgroup{iG};
        opt(indBG).iG=iG;
        opt(indBG).BestFPNch=showBestFPN;
        opt(indBG).PE=median(showBestFPN,1);
        opt(indBG).var=var(showBestFPN,1)/2;
        opt(indBG).tak=showBestFPN;
        opt(indBG).mean=mean(showBestFPN,1);
        opt(indBG).hp=mln_findpk(showBestFPN,1,Nbins);
        savefilename=[basicdir,'opt_para_FRN_hidden_HH_S',num2str(nc),'.mat'];
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