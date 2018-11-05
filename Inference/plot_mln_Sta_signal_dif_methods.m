function plot_mln_Sta_signal_dif_methods()

%dir='/Volumes/HWMac1/MULAN2/SPM/vnERP30GCBTED/ToutResults/';
dir='/Users/huifangwang/MULANII/programs/';
prenom='N30L10nmmCS100S1N4100';
datafile=[dir,'Sta_',prenom];
load(datafile,'mMULAN','xinfo','iMULAN');
load([dir,'Tout_',prenom],'Connectivity');

Theta=0.9:0.02:0.99;
Nth=length(Theta);
%xinfo.BM
dtoxaxis={'G12','G3','G4'};
index_iMULAN={'Chan','foo','G12','G3','G4','Nboot'};
[indxaxis,AUCtoshowind]=mln_chs_ind(index_iMULAN,dtoxaxis);
sizeMULAN=size(iMULAN.Mat);
subtoxaxis=sizeMULAN(indxaxis);

NGroup=prod(subtoxaxis);
Methsim=NaN(NGroup+1,Nth);

for iG=1:Nth
    [i12,i3,i4]=ind2sub(subtoxaxis,2);
    iMat=squeeze(iMULAN.Mat(:,:,i12,i3,i4,:));
for ith=1:Nth
[mlnMat,sim]=plot_mln_Sta_finalNetwork_im(iMat,Connectivity,Theta(ith));
Methsim(iG,ith)=sim;
end
end

for ith=1:Nth
[mlnMat,sim]=plot_mln_Sta_finalNetwork_im(mMULAN.Mat,Connectivity,Theta(ith));
Methsim(NGroup+1,ith)=sim;
end

Methsim
figure;
imagesc(Methsim);
