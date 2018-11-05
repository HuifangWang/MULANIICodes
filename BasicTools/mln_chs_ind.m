function [ind,updateind]=mln_chs_ind(indwhole,chind)

Nchind=length(chind);
updateind=indwhole;
ind=NaN(Nchind,1);
for ichind=1:Nchind
ind(ichind)=find(strcmp(chind(ichind),indwhole));
updateind{ind(ichind)}=[];
end
updateind=updateind(~cellfun('isempty',updateind));

