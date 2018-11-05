% huifang Wang,
% May 30, 2014
function nM=mln_mean(M,ind)
Mlength=length(ind);
nM=M;
for id=1:Mlength
    nM=mean(nM,ind(id));
end
nM=squeeze(nM);