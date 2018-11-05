% Huifang Wang, 29/10/13

function [MatNet,istep]=mln_algorithm_mlnvcs(MisInput,threshold,kforMF,Fthreshold,methodlog)
%load('Mis_nmm1.mat','MisInput');
data=MisInput;
numMFs=2;
inmftype='gbellmf';
outmftype='constant';
%threshold=[0.5,0.5,0.75,0.75];
%fis=mln_initialfis(data, numMFs, inmftype, outmftype,threshold,kforMF);
fis=mln_initialfis(data,methodlog,numMFs, inmftype, outmftype,threshold,kforMF);
[Nlink,NmethodsInput]=size(MisInput);
if NmethodsInput==length(methodlog)
    Input_AL=[MisInput,zeros(Nlink,2)];
else
Input_AL=MisInput(:,1:end-1);
end
in_output=zeros(Nlink,1);

%Fthreshold=0.5;
Nnode=(4*Nlink + 1)^(1/2)/2 + 1/2;
oldMat=zeros(Nnode,Nnode);
%oldMat=mln_Nlink2Stru(MisInput(:,end));
istep=1;
allMat1=cell(1,40);
allMatNet=cell(1,40);
while istep<40
    %% for each link
    in_output=zeros(Nlink,1);
    for ilink=1:Nlink
        input=Input_AL(ilink,:);
        %output_data=MisInput(end);
        in_output(ilink)=mln_evalfis(input, fis);
%        disp(['input=', num2str(input),';', num2str(in_output(ilink)),';', num2str(RMat(ilink))])
    end
    MatNet=mln_Nlink2Stru(in_output);
    iMat=MatNet.*(MatNet>=Fthreshold);
    inds=mln_sam_criteria(allMat1,istep,iMat);
    if isempty(inds)
        allMat1{istep}=iMat;
        allMatNet{istep}=MatNet;
        Input_AL=mln_updateInput(Input_AL,iMat);
        istep=istep+1;
    elseif inds==istep-1
            break;
    else
        sameMat=allMat1(inds:istep-1);
        iMat=mln_mean_cell(sameMat);
        MatNet=mln_mean_cell(allMatNet(inds:istep-1));
        if mln_Sim_Cos(iMat,allMat1{istep-1})>1-eps
            break;
        else
        allMat1{istep}=iMat;
        allMatNet{istep}=MatNet;
        Input_AL=mln_updateInput(Input_AL,iMat);
        istep=istep+1;
        end
    end
end
%save('Matfile.mat','allMat1');

function inds=mln_sam_criteria(allMat,istep,iMat)
inds=[];
if istep>1
    for im=1:istep-1
        if mln_Miscriteria(allMat{im},iMat)
            inds=im;
            break;
        end
    end
end
