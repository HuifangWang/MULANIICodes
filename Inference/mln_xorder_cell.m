%% to put the index of nCells inside of oCells
% Huifang Wang, Dec 14, 2013
function indx=mln_xorder_cell(oCells,nCells)
Nindx=length(nCells);
indx=NaN(Nindx,1);
for iindx=1:Nindx
indx(iindx)=find(strcmp(nCells{iindx},oCells));
end