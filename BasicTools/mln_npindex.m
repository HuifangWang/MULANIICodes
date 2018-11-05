function indbest=mln_npindex(vcs,bestvcs)
lenvcs=length(bestvcs);
indbest=NaN(lenvcs,1);
for ivcs=1:lenvcs
    indbest(ivcs,1)=find(vcs>=bestvcs(ivcs),1);
end