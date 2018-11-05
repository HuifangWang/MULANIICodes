function [isim,icos]=mln_sim_ngbmat(iMat,ngb_Mat)
% this function is used to get the everage similarity value and cos similarity
% value from iMat to it neighbor Mats

Ngb=length(ngb_Mat);
isimV=NaN(1,Ngb);
iCosV=NaN(1,Ngb);
for igb=1:Ngb
    isimV(1,igb)=mln_similiarity2M(iMat,ngb_Mat{igb});
    iCosV(1,igb)=mln_Sim_Cos(iMat,ngb_Mat{igb});
end
isim=mean(isimV);
icos=mean(iCosV);
