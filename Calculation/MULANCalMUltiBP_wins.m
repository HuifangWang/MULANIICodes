function MULANCalMUltiBP_wins(dirname,prenom,paramsfile,wins,VGroupMethlog)
% Huifang Wang, Marseille, August 26, 2013, Calculate all datasets with
% basic prenom
if ischar(paramsfile)
paramsfile=[dirname,'/',paramsfile];
load(paramsfile,'calParams');
else
    calParams=paramsfile;
end

calParams.defwindow=wins;
MulanCal(dirname,prenom,calParams,VGroupMethlog);
mln_Result2file(dirname,prenom,VGroupMethlog)