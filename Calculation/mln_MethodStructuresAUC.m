% Huifang Wang, Sep. 26, 2013
% Huifang Wang, Nov. 26, 2013 update Method Structure AUC for mlnI
% Mar 31, 2014 added AUC
% September 22 for Cos Similiarity
function mln_MethodStructuresAUC(dirname,prenom)

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
    %% for MULAN results for different thresholds
    SimTheta=0.75:0.02:0.96;
    Nstheta=length(SimTheta);
    neededfilemln=['./',dirname,'/ToutResults/mlnRes_',prenom,'.mat'];
    if exist(neededfilemln,'file')
        load(neededfilemln,'iMULAN','mMULAN','param');
        G12=param.methodlog.G12;
        G3=param.methodlog.G3;
        G4=param.methodlog.G4;
        th1=param.Thres.th1;
        th2=param.Thres.th2;
        mlnth=param.Thres.mth;
        
        NiG1=size(G12,1);
        NiG3=length(G3);
        NiG4=length(G4);
        Nith1=length(th1);
        Nith2=length(th2);
        Nmlnth=length(mlnth);
        MatRef=param.Rmat;
        iMULAN.AUC=NaN(NiG1,NiG3,NiG4,Nith1,Nith2,Nmlnth);
        iMULAN.MSim=NaN(NiG1,NiG3,NiG4,Nith1,Nith2,Nmlnth,Nstheta);
        iMULAN.COS=NaN(NiG1,NiG3,NiG4,Nith1,Nith2,Nmlnth,Nstheta);
        Methods='MULAN';
        %iMULAN.MSim=iMULAN.AUC;
        for iG1=1:NiG1
            for iG3=1:NiG3
                for iG4=1:NiG4
                    for ith1=1:Nith1
                        for ith2=1:Nith2
                            for imlnth=1:Nmlnth
                                
                                iMat=iMULAN.Mat(:,:,iG1,iG3,iG4,ith1,ith2,imlnth);
                                
                                
                                for isth=1:Nstheta
                                    iMat=iMat.*(iMat>SimTheta(isth));
                                    
                                    iMULAN.MSim(iG1,iG3,iG4,ith1,ith2,imlnth,isth)=mln_similiarity2M(iMat,MatRef);
                                    iMULAN.COS(iG1,iG3,iG4,ith1,ith2,imlnth,isth)=mln_Sim_Cos(MatRef,iMat);

                                    
                                end
                                [~,~,~,~,~,~,~,~,auc,~,~] =mln_calc_FalseRate(iMat,MatRef,mln_issymetricM(Methods),1);
                                iMULAN.AUC(iG1,iG3,iG4,ith1,ith2,imlnth)=auc;
                               % iMULAN.chis(iG1,iG3,iG4,ith1,ith2,imlnth)=auc;
                            end
                        end
                    end
                end
            end
        end
        
        %% for average
        for ith1=1:Nith1
            for ith2=1:Nith2
                for imlnth=1:Nmlnth
                    iMat=squeeze(mMULAN.Mmat(:,:,ith1,ith2,imlnth));
                    
                    [~,~,~,~,~,~,~,~,auc,~,~] =mln_calc_FalseRate(iMat,MatRef,0,1);
                    mMULAN.AUC(ith1,ith2,imlnth)=auc;
                    for isth=1:Nstheta
                        iMat=iMat.*(iMat>SimTheta(isth));
                        
                        mMULAN.MSim(ith1,ith2,imlnth,isth)=mln_similiarity2M(iMat,MatRef);
                        mMULAN.MSim(ith1,ith2,imlnth,isth)=mln_Sim_Cos(MatRef,iMat);
                        
                    end
                end
            end
            
            
        end
    end
    save(tosavefile,'mMULAN','iMULAN','param','-append');
end



