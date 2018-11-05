% huifang Wang, 10/29/13 this function is to get the statistically results
% for a set of given parameters
function outputname=mln_Stat_Valid_thv_mlnvcs(dir,prenom,threshold1,threshold2,kforMF,Fthreshold,nboot,nboot2,BM)
%Nchan=5;
%dir='/Volumes/HWMac1/MULAN2/SPM/vnERP30GCBTED/ToutResults/';
%dir='/Users/huifangwang/MULANII/programs/';
%prenom='N30L10nmmCS100S1N4100';

dir=[dir,'/ToutResults/'];
if nargin<3
    threshold1=0.6;
    threshold2=0.7;

else
    if ischar(threshold1)
        threshold1=str2double(threshold1);
        threshold2=str2double(threshold2);
    end
   
end
 threshold=[threshold1,threshold1,threshold2,threshold2];

if nargin<5
   kforMF=5;
else if ischar(kforMF)
        kforMF=str2double(kforMF);    
    end
end

if nargin<6
Fthreshold=0.7; % thm in the paper
else if ischar(Fthreshold)
        Fthreshold=str2double(Fthreshold);    
    end
end

if nargin<7
nboot=1000;
end
if nargin<8
nboot2=400;
end
dataprenom=[dir,'Tout_',prenom];
if nargin<9 
BM.G12={'BCorrU','PCorrU';'BTEU','PTEU'}; % I tried at for S1-S5
%BM.G12={'BH2U','PH2U';'BTEU','PTEU'};  % tired it on S6-S10
BM.G3={'BTED'};
BM.G4={'hmvar','ffDTF'};
end
G12=BM.G12;
G3=BM.G3;
G4=BM.G4;
load(dataprenom,'Connectivity');
if exist('Connectivity','var')
xinfo.MatRef=Connectivity;
end
xinfo.BM=BM;
xinfo.threshold=threshold;
xinfo.kforMF=kforMF;
xinfo.Fthreshold=Fthreshold;
xinfo.nboot=nboot;
xinfo.nboot2=nboot2;

th1name=num2str(threshold1*100);
th2name=num2str(threshold2*100);
outputname=['Sta_',prenom,'th12',th1name,th2name,'.mat'];
savestatisfile=[dir,outputname];
if ~exist(savestatisfile,'file')
    
methods_used={G12{:},G3{:},G4{:}};
NiG1=size(G12,1);
NiG3=length(G3);
NiG4=length(G4);
%[InputMM,InputBS,RMat]=mln_generate_BS_Input(dataprenom,List4Methods,nboot);
InputMM=mln_GenerateTDM(dataprenom,methods_used,nboot,nboot2);
indG=[3,4,5];
for iboot=1:nboot2
    for iG1=1:NiG1
        for iG3=1:NiG3
            for iG4=1:NiG4
                methodlog=[G12(iG1,:),G3{iG3},G4{iG4}];
                
                %Ninput=length(methodlog);
                pause(0.01)
                indMM=InputMM.methods;
                %ind_methods=[1 3 5 6]';%
                ind_methods=mln_chs_ind(indMM,methodlog);
                iInput=InputMM.values(:,ind_methods,iboot);
                [iMat,istep]=mln_algorithm_mlnvcs(iInput,threshold,kforMF,Fthreshold,methodlog);
                iMULAN.Mat(:,:,iG1,iG3,iG4,iboot)=iMat;
                iMULAN.step(:,:,iG1,iG3,iG4,iboot)=istep;
            end
        end
    end
    allMatiboot=iMULAN.Mat(:,:,:,:,:,iboot);
    
    mMULAN.Mat(:,:,iboot)=mln_mean(allMatiboot,indG);
end


save(savestatisfile,'iMULAN','mMULAN','xinfo');
end

function Input=mln_GenerateTDM(dataprenom,methods_used,nboot,nboot2)
MatA=load(dataprenom,methods_used{:});
Nchan=size(MatA.(methods_used{1}),1);
Nlink=Nchan*(Nchan-1);
Nmethod=length(methods_used);
Input.values=NaN(Nlink,Nmethod,nboot2);
Input.methods=methods_used;

nbins=100;
Norder=2;

for imethod=1:Nmethod
    Mat1=MatA.(methods_used{imethod});
    
    [Nnode,~,Nwin]=size(Mat1);
    isizeboot=Nlink;
    
    MatBoot=NaN( isizeboot*nboot,1);
    
    iMat=abs(mean(Mat1,3));
    iMat(1:size(iMat,1)+1:end) = NaN;
    a=reshape(iMat,Nnode*Nnode,1);
    MatBoot(1:isizeboot,1)=a(~isnan(a));
    
    for iboot=2:nboot
        iBoot_ind=randi(Nwin,[1,Nwin]);
        iBoot_Mat=Mat1(:,:,iBoot_ind);
        iMat=abs(mean(iBoot_Mat,3));
        iMat(1:size(iMat,1)+1:end) = NaN;
        a=reshape(iMat,Nnode*Nnode,1);
        MatBoot((iboot-1)*isizeboot+1:iboot*isizeboot,1)=a(~isnan(a));
    end
    xi=linspace(min(MatBoot),max(MatBoot),nbins);
    n_elements = histc(MatBoot,xi);
    c_elements = cumsum(n_elements);
    c_elements=c_elements./max(c_elements);
    %% seperate MatBoot to MBLink=M(Nlink,)
    %% the means for all windows
    MBLink=reshape(MatBoot,Nlink,nboot);
    for iboot=1:nboot2
        for ipairs=1:Nlink
            Input.values(ipairs,imethod,iboot)=c_elements(find(xi>=MBLink(ipairs,iboot),1));
        end
    end
    %         %% the means for  N-1 windows
    %         Vjack=1:Nwin;
    %         for ijackknife=1:Nwin
    %             a=Vjack;
    %             a(ijackknife)=[];
    %
    %         lMat=abs(mean(Mat1(:,:,a),3));
    %         lMat(1:size(iMat,1)+1:end) = NaN;
    %         a=reshape(lMat,Nnode*Nnode,1);
    %         pairs=a(~isnan(a));
    %         for ipairs=1:Nlink
    %         iTDM(ipairs+ijackknife*Nlink,imethod)=c_elements(find(xi>=pairs(ipairs),1));
    %         end
    %
    %         end
    %
end

function List4Methods=mln_GroupMethods(G1,G3,G4)
List4Methods={};
for iG1=1:length(G1)
    for iG3=1:length(G3)
        for iG4=1:length(G4)
            if isempty(List4Methods)
                List4Methods{1}=[G1(iG1,:),G3{iG3},G4{iG4}];
            else
                List4Methods{end+1}=[G1(iG1,:),G3{iG3},G4{iG4}];
            end
        end
    end
end
