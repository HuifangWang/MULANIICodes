function mln_save_para_gb_vcs_hidden_hc_allinall()

%% save the calculated parameters
%vnc=[20,30];
nc=84;
ncrit=0;

Nstas=50;
NP=5100;
methodgroup={{'BTEU','BCorrU','hmvar','ffDTF'},{'BCorrU','BTEU','BTED','ffDTF'},...,
    {'BCorrU','PCorrU','hmvar','ffDTF'},{'BCorrU','PCorrU','BTED','hmvar'},{'BCorrU','COH1','BCorrD','PDC'}};
BG=[1,2,3,4,5];
%BG=[2,5];
%basicdir='/Users/huifangwang/MULANa/humanCon/HCRN84S100N/hidden/Win1500/ToutResults/';
      
basicdir='/Users/huifangwang/MULANa/humanCon/HC84S100hiT/hidden/';
filesubdir='Win1500/ToutResults/';
Nbins=20;

bgv=nchoosek(1:6,2);
%
ds=0;
nfile=0;

        iwins=1500;
        
        
    for indBG =1:length(BG)
         iG=BG(indBG);
         showBestFPN=[];
for ibgv=1:length(bgv)
    

       
        
    
        
        
        fitness='FPN';
        for istru =1:50 %istru =[1:4,10:49]%1:Nstas
            datasets=['erpN',num2str(nc),'S',num2str(istru),'N',num2str(NP),'b',num2str(bgv(ibgv,1)),num2str(bgv(ibgv,2))];
            filename=[basicdir,filesubdir,'Fit_',fitness,'G',num2str(iG),'_',num2str(iwins),'_',datasets];
            if ds>0
                filename=[filename,'ds',num2str(ds)];
            end
            
            
                try
                
                    
                load(filename,'state');
                catch
                    continue
                end
                nfile=1+nfile;
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
        savefilename=[basicdir,'opt_para_FRN_hidden_allin1',num2str(nc),'.mat'];
        save(savefilename,'opt')
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