function [aucws,thetaws,FPws,FNws,FPNws,bithws]=mln_AUC_TPN_vcs(iMat,gtMat,wscs)

%% this AUC based on sample and directed matrix

aucws=NaN(1,2);
thetaws=NaN(1,2);
FPws=NaN(1,2);
FNws=NaN(1,2);
FPNws=NaN(1,2);
bithws=NaN(1,2);

for iws=1:2
bigtMat=gtMat(:,:)>=wscs(iws);

H0=iMat(~bigtMat);
H1=iMat(bigtMat);
NH0=length(H0);
NH1=length(H1);

%% 
FPv=NaN(1,3);
FNv=NaN(1,3);
thetav=[min(H1),max(H0),mean([min(H1),max(H0)])];
for ith=1:3
FPv(ith)=length(H0(H0>=thetav(ith)));
FNv(ith)=length(H1(H1<thetav(ith)));
end

[FPNws(iws),bith]=min(FPv+FNv,[],2);
FPws(iws)=FPv(bith);
FNws(iws)=FNv(bith);
thetaws(iws)=thetav(bith);
bithws(iws)=bith;
%% AUC
Nsteps=100;
ih=linspace(min(min(iMat)),max(max(iMat)),Nsteps);
t=ih(ih>=0);
Nt=length(t);

Fpr=zeros(Nt,1);
Tpr=zeros(Nt,1);
Fpr(1)=1;
Tpr(1)=1;
for i=2:Nt
    Fpr(i)=(sum(H0>=t(i)))/NH0;
    Tpr(i)=(sum(H1>=t(i)))/NH1;
end

if isempty(find(iMat))
    aucws(iws)=0.5;
else

auc = 0.5*sum( (Fpr(2:end)-Fpr(1:end-1)).*(Tpr(2:end)+Tpr(1:end-1)) );
aucws(iws) = abs(auc);
end

end

