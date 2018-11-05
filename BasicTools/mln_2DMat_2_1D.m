function MatBoot=mln_2DMat_2_1D(iMat)
Nnode=size(iMat,1);
iMat(1:Nnode+1:end) = NaN;
    a=reshape(iMat,Nnode*Nnode,1);
    MatBoot=a(~isnan(a));