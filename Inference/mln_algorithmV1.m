function [iMat,istep]=mln_algorithmV1(MisInput,threshold,kforMF)
%load('Mis_nmm1.mat','MisInput');
data=MisInput;
numMFs=2;
inmftype='gbellmf';
outmftype='constant';
%threshold=[0.5,0.5,0.75,0.75];
fis=mln_initialfis(data, numMFs, inmftype, outmftype,threshold,kforMF);
%fis=initialMLNfisV1(data, numMFs, inmftype, outmftype,threshold,kforMF);
Nlink=size(MisInput,1);
Input_AL=MisInput(:,1:end-1);
in_output=zeros(Nlink,1);
threshold=0.5;
Nnode=(4*Nlink + 1)^(1/2)/2 + 1/2;
oldMat=zeros(Nnode,Nnode);
%oldMat=mln_Nlink2Stru(MisInput(:,end));
istep=1;
while istep<40
%% for each link
for ilink=1:Nlink
input=Input_AL(ilink,:);
%output_data=MisInput(end);
in_output(ilink)=mln_evalfis(input, fis);
end
MatNet=mln_Nlink2Stru(in_output);
iMat=MatNet.*(MatNet>=threshold);

if ~mln_Miscriteria(oldMat,iMat)
Input_AL=mln_updateInput(Input_AL,iMat);
oldMat=iMat;
istep=istep+1;
else
    break;
end
end
%save('errorfile.mat','errOutput');
