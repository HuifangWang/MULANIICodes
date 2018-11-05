% huifang Wang, Jan, 30, 2015 This function outputs mln algorithm
% connection matrix for given nboot2, BMG, thresholds12
function outputname=mln_alg_output(dir,prenom,threshold1,threshold2,kforMF,Fthreshold,nboot,nboot2,BMG,add)


dir=[dir,'/ToutResults/'];
if nargin<3
    threshold1=0.5:0.05:0.6;
    threshold2=0.6:0.05:0.7;
    
else
    if ischar(threshold1)
        threshold1=str2double(threshold1);
        threshold2=str2double(threshold2);
    end
    
end

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
    nboot2=20;
end
dataprenom=[dir,'Tout_',prenom];
if nargin<9
    BMG={{'BCorrU','PCorrU','BCorrD','hmvar'},{'BTEU','PTEU','BCorrD','ffDTF'}};
end
if nargin<10
    add='';
else
add=['_',add];
end

    mmRule=[2,0.9,0.9;1,0.8,0.8];

outputname=['Sta_',prenom,add,'.mat'];
savestatisfile=[dir,outputname];
if ~exist(savestatisfile,'file')
    
    
    xinfo.th1=threshold1;
    xinfo.th2=threshold2;
    xinfo.kforMF=kforMF;
    xinfo.Fthreshold=Fthreshold;
    xinfo.nboot=nboot;
    xinfo.nboot2=nboot2;
    
    xinfo.BMG=BMG;
    NmethodG=length(BMG);
    Nth1=length(threshold1);
    Nth2=length(threshold2);
    %iMULAN.Mat=NaN(Node,Node,Nth1,Nth2,NmethodG,nboot2);
    for ith1=1:Nth1
        for ith2=1:Nth2
            
            threshold=[threshold1(ith1),threshold1(ith1),threshold2(ith2),threshold2(ith2)];
            
            for img=1:NmethodG
                methodlog=BMG{img};
                InputMM=mln_GenerateTDM(dataprenom,methodlog,nboot,nboot2);
                %indG=[3,4,5];
                
                for iboot=1:nboot2
                    
                    %Ninput=length(methodlog);
                    pause(0.01)
                    indMM=InputMM.methods;
                    %ind_methods=[1 3 5 6]';%
                    ind_methods=mln_chs_ind(indMM,methodlog);
                    iInput=InputMM.values(:,ind_methods,iboot);
                    [iMat,istep]=mln_algorithm_mlnvcs(iInput,threshold,kforMF,Fthreshold,methodlog);
                    iMULAN.Mat(:,:,ith1,ith2,img,iboot)=iMat;
                    iMULAN.step(:,:,ith1,ith2,img,iboot)=istep;
                end
            end
            
        end

    end
    Nchan=size(iMULAN.Mat,1);
    allGMat=reshape(iMULAN.Mat,[Nchan,Nchan,Nth1*Nth2*NmethodG*nboot2]);
    [MULAN.mdMat,MULAN.mnMat,MULAN.fMat]=mln_output_sta_CM(allGMat,mmRule);
    load(dataprenom,'Connectivity');
    if exist('Connectivity','var')
        xinfo.MatRef=Connectivity;
        
        [Eva.AUC(1),Eva.th(1)]=mln_AUC(MULAN.mdMat,xinfo.MatRef,0,1);
        [Eva.AUC(2),Eva.th(2)]=mln_AUC(MULAN.mnMat,xinfo.MatRef,0,1);
        if Eva.AUC(1)>Eva.AUC(2)
            MULAN.fMat=MULAN.mdMat.*(MULAN.mdMat>Eva.th(1));
        else
            MULAN.fMat=MULAN.mnMat.*(MULAN.mnMat>Eva.th(2));
        end
        iMatRef=Connectivity.*(Connectivity>=0.7);%%!Attension    
        Eva.COS=mln_Sim_Cos(iMatRef,MULAN.fMat);
    save(savestatisfile,'iMULAN','MULAN','Eva','xinfo');
    else
        save(savestatisfile,'iMULAN','MULAN','xinfo');
    end
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

end






