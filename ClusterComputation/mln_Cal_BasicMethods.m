function mln_Cal_BasicMethods(dirname,dataname,paramsfile,VGroupMethlog)

% Huifang Wang, Marseille, Sep.5 2014, Calculate all datasets with
% basic prenom
if ischar(paramsfile)
paramsfile=[dirname,'/',paramsfile];
load(paramsfile,'calParams');
else
    calParams=paramsfile;
end

mln_Basic_Cal_Detail(dirname,dataname,calParams,VGroupMethlog);
mln_Result2file(dirname,dataname,VGroupMethlog)