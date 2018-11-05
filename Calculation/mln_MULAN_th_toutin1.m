%% Huifang Wang, Marseille, 02/11/13

function mln_MULAN_th_toutin1(dirname,prenom,Methodthresholdfile)
% Example: mlnC_MULAN_th 5 nmm TPnmmCS100S24N3000 SetMethodsTh

load([dirname,'/',Methodthresholdfile,'.mat'],'G12','G3','G4','thresholdsM');
resfilename=[dirname,'/ToutResults','/Tout_',prenom,'.mat'];
th1=thresholdsM.th1;
th2=thresholdsM.th2;
mlnth=thresholdsM.mlnth;
kforMF=5;
%Fthreshold=str2double(Fthreshold);
tosavefile=[dirname,'/ToutResults','/mlnRes_',prenom,'.mat'];
if exist(tosavefile,'file')
    [cdata] = imread('MULANLOGO.png');
msgbox('Calculation already there','Existed','custom',cdata);
else

NiG1=size(G12,1);
NiG3=length(G3);
NiG4=length(G4);
Nith1=length(th1);
Nith2=length(th2);
Nmlnth=length(mlnth);

samp_methods=G3{1};
samp=load(resfilename,samp_methods);
Nchan=size(samp.(samp_methods),1);
%iMULAN.AUC=NaN(NiG1,NiG3,NiG4,Nith1,Nith2,Nmlnth);
iMULAN.Mstep=NaN(NiG1,NiG3,NiG4,Nith1,Nith2,Nmlnth);
%iMULAN.MSim=iMULAN.AUC;
iMULAN.Mat=NaN(Nchan,Nchan,NiG1,NiG3,NiG4,Nith1,Nith2,Nmlnth);
savedsubfolder=[dirname,'/ToutResults/', prenom];



            for ith1=1:Nith1
                for ith2=1:Nith2
                    for imlnth=1:Nmlnth
                    
                    getonefile=[savedsubfolder,'/mlnRes_',prenom,'th1',num2str(th1(ith1)*100),'th2',num2str(th2(ith2)*100),'mth',num2str(mlnth(imlnth)*10),'.mat'];
                    
                    singleth=load(getonefile);
                    
                    iMULAN.Mat(:,:,:,:,:,ith1,ith2,imlnth)=singleth.iMULAN.Mat;
                    iMULAN.Mstep(:,:,:,ith1,ith2,imlnth)=singleth.iMULAN.Mstep;
                    %iMULAN.MSim(iG1,iG3,iG4,ith1,ith2,imlnth)=Sim;
                    
                    %iMULAN.AUC(iG1,iG3,iG4,ith1,ith2,imlnth)=auc;
                    end
                end
            end
  

%% Average Values
      mMULAN.Mmat=squeeze(mean(mean(mean(iMULAN.Mat,5),4),3));
iRMat=load(resfilename,'Connectivity');

param.methodlog.G12=G12;
param.methodlog.G3=G3;
param.methodlog.G4=G4;
if ~isempty(iRMat) && ~isempty(fields(iRMat))
param.Rmat=iRMat.Connectivity;
end
param.Thres.th1=th1;
param.Thres.th2=th2;
param.Thres.mth=mlnth;
%param.Thres.SimTheta=SimTheta;

save(tosavefile,'param','iMULAN','mMULAN');
    [cdata] = imread('MULANLOGO.png');
msgbox('Mulan calculation completed','Success','custom',cdata);
end



