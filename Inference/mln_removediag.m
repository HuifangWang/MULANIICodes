function MatC=mln_removediag(MatC)
Nnode=size(MatC,1);
MatC(1:size(MatC,1)+1:end) = NaN;
a=reshape(MatC,Nnode*Nnode,1);
MatC=a(~isnan(a));