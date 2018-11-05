% Huifang Wang, Sep. 26, 2013
% Huifang Wang, Nov. 26, 2013 update Method Structure AUC for mlnI
% Mar 31, 2014 added AUC
function mln_AUC_Tout(dirname,prenom)

neededfile=['./',dirname,'/ToutResults/Tout_',prenom,'.mat'];
tosavefile=['./',dirname,'/ToutResults/AUC_',prenom,'.mat'];
if ~exist(tosavefile,'file')
    
    if exist(neededfile,'file')
            
        calresult=load(neededfile);
        fieldname=fieldnames(calresult);
        % remove params and connectivity
        fieldname(strcmp('Params',fieldname))=[];
        
        fieldname(strcmp('Connectivity',fieldname))=[];
        %%
        %idx=find(strcmp('pCOH1',fieldname));
        %fieldname(idx)=[];
        %%
        Nmethod=length(fieldname);
        Meths.MSAUC=zeros(Nmethod,1);
        Meths.Mth=zeros(Nmethod,1);
        
        
        Meths.Methodnames=fieldname;
        
        Meths.Connectivity=calresult.Connectivity;
        
        for imethod=1:Nmethod
            Methods=fieldname{imethod};
            Mat=calresult.(Methods);
            if isnan(max(Mat))
                
                auc=0;
                chis=1;
                
            else
                iMat=mean(abs(Mat),3);
                iMat=iMat-diag(diag(iMat));
                [~,~,~,~,~,~,~,~,auc,~,thresh3] =mln_calc_FalseRate(iMat,calresult.Connectivity,mln_issymetricM(Methods),1);
                chis=mln_chis(Mat,thresh3(1));
            end
            Meths.MSAUC(imethod,1)=auc;
            
            if isempty(chis)
                pause
            end
            Meths.Mth(imethod,1)=chis;
            Meths.Mat(:,:,imethod)=iMat;
        end
        
        
        save(tosavefile,'Meths');
    end
end
if exist(tosavefile,'file')
    saveddata=load(tosavefile);
   if ismember(fieldnames(saveddata),'mMULAN')
       return
   end
end
   
end



