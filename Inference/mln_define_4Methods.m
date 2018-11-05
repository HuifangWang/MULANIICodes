%Huifang Nov. 13
function List4Methods=mln_define_4Methods(Methodfile)
basicMethod(1,:)={'BCohF','PCohF'}; %1
basicMethod(2,:)={'BCorrU','PCorrU'}; %2
basicMethod(3,:)={'BH2U','PH2U'}; %3
% basicMethod(4,:)={'BCorrD','PCorrD'};%4
% basicMethod(5,:)={'BH2D','PH2D'};%5
% basicMethod(6,:)={'COH1','pCOH1'};%6
% basicMethod(7,:)={'COH2','pCOH2'};%7
% basicMethod(8,:)={'BMITU','PMITU'};%8
% basicMethod(9,:)={'BMITD1','PMITD1'};
% basicMethod(10,:)={'BMITD2','PMITD2'};
% basicMethod(11,:)={'BCohW','PCohW'};

G3={'GC','CondGC'};

%G4={'oPDCF','GPDC','PDC','DTF','DC1','GPDC','dDTF','ffDTF','GGC'};
G4={'oPDCF','GPDC'};

List4Methods={};
for iG1=1:length(basicMethod)
    for iG3=1:length(G3)
        for iG4=1:length(G4)
if isempty(List4Methods)
            List4Methods{1}=[basicMethod(iG1,:),G3{iG3},G4{iG4}]; 
else
            List4Methods{end+1}=[basicMethod(iG1,:),G3{iG3},G4{iG4}]; 
end
        end
    end
end

%save([Methodfile,'SetMethods.mat'],'basicMethod','G3','G4','thresholdsM');
