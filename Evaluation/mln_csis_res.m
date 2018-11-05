function mln_csis_res()
%% this function is for get the consistence matrix from MULAN results files
dirname='/Volumes/HWMac1/MULAN2/SPM/vnERP';
prenom='N30L10nmmCS80S7N4100';

%% for find the consistence with vaying of thresholds
neededfilemln=[dirname,'/ToutResults/mlnRes_',prenom,'.mat'];
tosavefile=[dirname,'/ToutResults/cons_',prenom,'.mat'];
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
    iMULAN.Sim1=NaN(NiG1,NiG3,NiG4,Nith1,Nith2,Nmlnth); % Similarity
    iMULAN.Cos2=NaN(NiG1,NiG3,NiG4,Nith1,Nith2,Nmlnth); % Cosine
    
    %iMULAN.MSim=iMULAN.AUC;
    for iG1=1:NiG1
        for iG3=1:NiG3
            for iG4=1:NiG4
                for ith1=1:Nith1
                    for ith2=1:Nith2
                        for imlnth=1:Nmlnth
                            
                            iMat=iMULAN.Mat(:,:,iG1,iG3,iG4,ith1,ith2,imlnth);
                            
                            
                            % get neighbor Mat
                            if ith1==1
                                indx1=1:2;
                            else if ith1==Nith1
                                    indx1=ith1-1:Nith1;
                                else
                                    indx1=ith1-1:ith1+1;
                                end
                            end
                            
                            if ith2==1
                                indx2=1:2;
                            else if ith2==Nith2
                                    indx2=ith2-1:ith2;
                                else
                                    indx2=ith2-1:ith2+1;
                                end
                            end
                            nngk=0; % count of neighbors
                            for iindx1=1:length(indx1)
                                for iindx2=1:length(indx2)
                                    if ~(indx1(iindx1)==ith1 && indx2(iindx2)==ith2)
                                        ngb_Mat{nngk+1}=iMULAN.Mat(:,:,iG1,iG3,iG4,indx1(iindx1),indx2(iindx2),imlnth);
                                        nngk=nngk+1;
                                    end
                                end
                            end
                            [isim1,icos1]=mln_sim_ngbmat(iMat,ngb_Mat);
                            iMULAN.Sim1(iG1,iG3,iG4,ith1,ith2,imlnth)=isim1; % Similarity
                            iMULAN.Cos2(iG1,iG3,iG4,ith1,ith2,imlnth)=icos1;
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
                
                if ith1==1
                    indx1=2;
                else if ith1==Nith1
                        indx1=ith1-1;
                    else
                        indx1=ith1-1:ith1+1;
                    end
                end
                
                if ith2==1
                    indx2=2;
                else if ith2==Nith2
                        indx2=ith2-1;
                    else
                        indx2=ith2-1:ith2+1;
                    end
                end
                nngk=0; % count of neighbors
                for iindx1=1:length(indx1)
                    for iindx2=1:length(indx2)
                        if ~(indx1(iindx1)==ith1 && indx2(iindx2)==ith2)
                                        ngb_Mat{nngk+1}=mMULAN.Mmat(:,:,indx1(iindx1),indx2(iindx2),imlnth);
                                        nngk=nngk+1;
                        end
                    end
                end
                [isim1,icos1]=mln_sim_ngbmat(iMat,ngb_Mat);
                mMULAN.Sim1(ith1,ith2,imlnth)=isim1; % Similarity
                mMULAN.Cos2(ith1,ith2,imlnth)=icos1;
            end
        end
    end
    
end


      save(tosavefile,'mMULAN','iMULAN','param');
end



