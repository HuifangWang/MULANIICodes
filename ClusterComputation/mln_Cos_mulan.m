   function mln_Cos_mulan(dirname,prenom)
   % Calculate the Cos similarity of mulan
   if dirname(1)=='/'
       dirname=dirname;
   else
       dirname=['./',dirname];
   end
   
   tosavefile=[dirname,'/ToutResults/COS_',prenom];
   if ~exist(tosavefile,'file')
    neededfilemln=[dirname,'/ToutResults/',prenom];
    if exist(neededfilemln,'file')
        load(neededfilemln,'iMULAN','mMULAN','xinfo');
        G12=xinfo.BM.G12;
        G3=xinfo.BM.G3;
        G4=xinfo.BM.G4;
        Theta=0.99;
        NiG1=size(G12,1);
        NiG3=length(G3);
        NiG4=length(G4);
        Methods='MULAN';
        incof=0.7:0.01:0.97;
        Nncof=length(incof);
        
        Nnode=size(iMULAN.Mat,1);
        MatRef=xinfo.MatRef;
        
        iMULAN.AUC=NaN(NiG1,NiG3,NiG4);
        iMULAN.MSim=NaN(NiG1,NiG3,NiG4,Nncof);
        iMULAN.mlnMat=NaN(Nnode,Nnode,NiG1,NiG3,NiG4,Nncof);
        %iMULAN.MSim=iMULAN.AUC;
        for iG1=1:NiG1
            for iG3=1:NiG3
                for iG4=1:NiG4
                    iMatiG=squeeze(iMULAN.Mat(:,:,iG1,iG3,iG4,:));
                    miMat=mean(iMatiG,3);
                    [~,~,~,~,~,~,~,~,auc,~,~] =mln_calc_FalseRate(miMat,MatRef,mln_issymetricM(Methods),1);
                    iMULAN.AUC(iG1,iG3,iG4)=auc;
                    % iMULAN.chis(iG1,iG3,iG4,ith1,ith2,imlnth)=auc;
                    for i_incof=1:Nncof
                        ncof=incof(i_incof);
                        [mlnMat,sim]=plot_mln_Sta_finalNetwork_Mat(iMatiG,MatRef,Theta,ncof);
                        iMULAN.MSim(iG1,iG3,iG4,i_incof)=sim;
                        iMULAN.mlnMat(:,:,iG1,iG3,iG4,i_incof)=mlnMat;
                        
                    end
                end
                
            end
        end
        
        
        
        
        %% for average
       
                    iMat=squeeze(mean(mMULAN.Mat,3));
                    
                    [~,~,~,~,~,~,~,~,auc,~,~] =mln_calc_FalseRate(iMat,MatRef,0,1);
                    mMULAN.AUC=auc;
                               for i_incof=1:Nncof
                                   ncof=incof(i_incof);
                    [mlnMat,sim]=plot_mln_Sta_finalNetwork_Mat(iMat,MatRef,Theta,ncof);
                               mMULAN.MSim(i_incof)=sim;
                               mMULAN.mlnMat(:,:,i_incof)=mlnMat;
                               end
                    save(tosavefile,'mMULAN','iMULAN','xinfo');          
                end
    
    
   end
   
   
   
   
function [mlnMat,sim]=plot_mln_Sta_finalNetwork_Mat(miMat,MatRef,Theta,ncof)

[Nnode,~,Nboots]=size(miMat);
xvalues=linspace(0,1,50);
xind=mln_v_find(xvalues,Theta,1);
mlnMat=NaN(Nnode);
for inode=1:Nnode
    for jnode=1:Nnode
        
        pij=squeeze(miMat(inode,jnode,:));
        [ph,xh]=hist(pij,xvalues);
         if mln_sumhist_more2(ph,xind,Nboots,ncof,'AND')
             mlnMat(inode,jnode)=1;
         else
             mlnMat(inode,jnode)=0;
         end
    end
end

sim=mln_Sim_Cos(MatRef,mlnMat);

function xind=mln_v_find(xvalues,Theta,nN)
Nth=length(Theta);
xind=NaN(Nth,nN);
for ith=1:Nth
    xind(ith,:)=find(xvalues>Theta(ith),nN);
end

function TorF=mln_sumhist_more2(ph,xind,Nboots,ncof,flaglogic)
switch flaglogic
    case 'OR'
        TorF=0;
    case 'AND'
        TorF=1;
end
Nth=length(xind);
for ith=1:Nth
    if sum(ph(xind(ith):end))>Nboots*ncof(ith)
        switch flaglogic
            case 'OR'
                TorF=1;
                return;
        end
    else
        switch flaglogic
            case 'AND'
                TorF=0;
                return;
        end
    end
    
end