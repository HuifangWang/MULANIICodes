%% Function to calculate the pairs of nodes which have the common source 
% at most common source =2 at the current version
% Huifang Wang, Aug. 29 2013, Marseille
function C1=mln_CommonSource(A,Norder)

Nchan=size(A,1);
C1=zeros(Nchan);
% build the tree
NT=0;
for icol=1:Nchan
idcol=find(A(:,icol));
if length(idcol)>1
    NT=NT+1;
    NodeTress(NT).L0=icol;
    NodeTress(NT).L1=idcol;
    C1(idcol,idcol)=1;
    if Norder>1
        for iidcol=1:length(idcol)
            iL2node=find(A(:,idcol(iidcol)));
            aa= find(~ismember(iL2node,NodeTress(NT).L0) &  ~ismember(iL2node,NodeTress(NT).L1));
            if ~isempty(aa)
                NodeTress(NT).L2(iidcol).brand=[idcol(idcol~=idcol(iidcol))',iL2node(aa)'];
                C1(NodeTress(NT).L2(iidcol).brand,NodeTress(NT).L2(iidcol).brand)=1;
            end
        end
    end
end
end

    
C1=C1-diag(diag(C1));


