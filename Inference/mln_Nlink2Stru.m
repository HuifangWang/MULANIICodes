function MatNet=mln_Nlink2Stru(CS)
Nlink=length(CS);
Nnode=(4*Nlink + 1)^(1/2)/2 + 1/2;
updateCS=zeros(Nlink+Nnode,1);
indexl=1:Nnode*Nnode;
indexdiag=1:Nnode+1:Nnode*Nnode;
indexl(indexdiag)=[];
try updateCS(indexl)=CS;
catch
    error('index')
end
MatNet=reshape(updateCS,[Nnode,Nnode]);

