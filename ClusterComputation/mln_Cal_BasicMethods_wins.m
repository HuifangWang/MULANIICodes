function mln_Cal_BasicMethods_wins(dirname,dataname,paramsfile,VGroupMethlog,add)

% Huifang Wang, Marseille, Sep.5 2014, Calculate all datasets with
% basic prenom
if ischar(paramsfile)
paramsfile=[dirname,'/',paramsfile];
load(paramsfile,'calParams');
else
    calParams=paramsfile;
end

mln_Basic_Cal_Detail(dirname,dataname,calParams,VGroupMethlog,add);
mln_Result2file_add(dirname,dataname,VGroupMethlog,add)