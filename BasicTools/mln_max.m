% huifang Wang,
% May 30, 2014
function nM=mln_max(M,ind)
Mlength=length(ind);
nM=M;
for id=1:Mlength
    nM=max(nM,[],ind(id));
end
nM=squeeze(nM);

