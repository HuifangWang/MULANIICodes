% Huifang Wang, Feb 24, 2015
% AUC for dynamics
function mln_AUC_Tout_dyn(dirname,prenom)

if dirname(1)~='/'
    dirname=['./',dirname];
end

neededfile=[dirname,'/ToutResults/Tout_',prenom,'.mat'];
tosavefile=[dirname,'/ToutResults/AUC_',prenom,'.mat'];
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
        Nwins=size(calresult.(fieldname{1}),3);
        Connectivity=calresult.Connectivity;
        Ndy=length(Connectivity.switchTimes);
        Meths.MSAUC=zeros(Nmethod,Ndy);
        Meths.Mth=zeros(Nmethod,Ndy);
        
        
        Meths.Methodnames=fieldname;
        
        Meths.Connectivity=calresult.Connectivity;
        
        tstfwins=mln_tstf_2wins(Connectivity.switchTimes,calresult.Params);
        tstfwins(Ndy,2)=Nwins;
        for imethod=1:Nmethod
            Methods=fieldname{imethod};
            Mat=calresult.(Methods);
            for idwin=1:Ndy
            if isnan(max(Mat))
                
                auc=0;
                %chis=1;
                
            else
                iMat=mean(abs(Mat(:,:,tstfwins(idwin,:))),3);
                iMat=iMat-diag(diag(iMat));
                auc=mln_AUC(iMat,Connectivity.strus(:,:,idwin),mln_issymetricM(Methods),1);
                %chis=mln_chis(Mat,thresh3(1));
            end
            Meths.MSAUC(imethod,idwin)=auc;
            
            %if isempty(chis)
                %pause
            %end
            %Meths.Mth(imethod,iwin)=chis;
            Meths.Mat(:,:,imethod,idwin)=iMat;
            end
        end
        
        
        save(tosavefile,'Meths');
        
    end

    
end
   
end



