function [thf,ipfit]=mln_updata_thf(iMat,fMat,rangethf,auc)

epsmy=0.01;
%
trithf = linspace(rangethf(1),rangethf(2),5);
trifit = mln_trifit(iMat,fMat,trithf);
[trifits,Indtri]=sort(trifit);
icount=1;
while abs(trifits(5)-trifits(4))>epsmy && abs(trifits(5)-auc)>epsmy && icount>10
    trithf = linspace(trithf(Indtri(4)),trithf(Indtri(5)),5);
    
    trifit = mln_trifit(iMat,fMat,trithf);
    [trifits,Indtri]=sort(trifit);
    icount=icount+1;
end
thf=trithf(Indtri(5));
ipfit=trifits(5);


function trifit = mln_trifit(iMat,fMat,trithf)
trifit = zeros(5,1);
for itri = 1: 5
    iMatRef=iMat.*(iMat>=trithf(itri));
    trifit(itri) =mln_Sim_Cos(iMatRef,fMat);
end
