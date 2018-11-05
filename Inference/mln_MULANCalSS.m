%function [iMat,MatRef,istep,Sim,auc]=mln_MULANCalSS(datafile,threshold,methodlog,kforMF,Fthreshold)
function [iMat,MatRef,istep]=mln_MULANCalSS(datafile,threshold,methodlog,kforMF,Fthreshold)
%Nchan=5;
%datafile='Tout_N30nmmCS80S2N4100';
% idata=13;
%threshold=[0.7,0.7,0.8,0.8];
%methodlog={'BCorrU','PCorrU','GC','hmvar'};
%kforMF=5;
%Fthreshold=0.5;
%methodlog={'BCorrD','PCorrD','GC','oPDCF'};
%methodlog={'BCorrD','PCorrD','GC','GPDC'};
%methodlog={'BCorrD','PCorrD','GC','PGC'};
%threshold=mln_getthreshold4Sim(methodlog,savedSfile);
%methodlog={'COH1','pCOH1','GC','oPDCF'};
%methodlog={'COH1','pCOH1','GC','oPDCF'};
% methodlog={'BCohF','PCohF','GC','oPDCF'};
% methodlog={'BCohF','PCohF','GC','GPDC'};
% methodlog={'BCohF','PCohF','GC','GGC'};
% methodlog={'BCohF','PCohF','GC','PGC'};
%dataprenom=[prenom,num2str(idata)];
MisInput=mln_generateinput4neurofuzzy(datafile,methodlog);
Connectivity=MisInput(:,end);
if ~isnan(Connectivity)
MatRef=mln_Nlink2Stru(Connectivity);
else
    MatRef=[];
end
 %threshold=[0.5,0.5,0.7,0.7];
       
    [iMat,istep]=mln_algorithm(MisInput,threshold,kforMF,Fthreshold);
    %Sim=mln_similiarity2M(iMat,MatRef);
 %[~,~,~,~,~,~,~,~,auc,~,~] =mln_calc_FalseRate(iMat,MatRef,0,1); 