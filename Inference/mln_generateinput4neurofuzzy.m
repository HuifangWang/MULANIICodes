%% function produce the training data for the neurofuzzy controller
%% Huifang Wang, Aug 24 2013, Marseille
%% Huifang Wang, Aug 28, 2013, Marseille, generate Nwin-1 more data
%% Huifang Wang, Aug 30, Marseille consider two order indirect connection
%% and common source

function MisInput=mln_generateinput4neurofuzzy(datafile,methodlog)

% this function will generate input for MULAN inference systems
%[x_1,x_2,x_3,x_4,C_1,C_2, Y], x_1 is the cdf of CS of each links, Y is the normalized reference link
%[0, 1],C_i is the two order indirect connection and common source (1/0),
%which is set to 0

% BCorrU;PCorrU;GC;oPDCF; Y is the reference


%dataprenom='nmm1';

%methodlog={'BMITU','PMITU','GC','oPDCF'};
nboot=1000;

MisInput=GenerateTDM(datafile,methodlog,nboot);


%savedfilename=['Mis_',dataprenom,'.mat'];
%save(savedfilename,'MisInput');

function iTDM=GenerateTDM(datafile,methodlog,nboot)
% read the files

nbins=100;
%datafile=['Tout_',dataprenom,'.mat'];
Ninput=length(methodlog);
MatA=load(datafile,methodlog{:});
[Nnode,~,Nwin]=size(MatA.(methodlog{1}));
Nlink=Nnode*(Nnode-1);
iTDM=zeros(Nlink,Ninput+3);

for imethod=1:Ninput
    clear a
    Mat1=MatA.(methodlog{imethod});
    
    [Nnode,~,Nwin]=size(Mat1);
    isizeboot=Nlink;
    
    MatBoot=NaN( isizeboot*nboot,1);
    
    iMat=mean(abs(Mat1),3);
    iMat(1:size(iMat,1)+1:end) = NaN;
    a=reshape(iMat,Nnode*Nnode,1);
    MatBoot(1:isizeboot,1)=a(~isnan(a));
    
    for iboot=2:nboot
        iBoot_ind=randi(Nwin,[1,Nwin]);
        iBoot_Mat=Mat1(:,:,iBoot_ind);
        iMat=mean(abs(iBoot_Mat),3);
        iMat(1:size(iMat,1)+1:end) = NaN;
        ia=reshape(iMat,Nnode*Nnode,1);
        MatBoot((iboot-1)*isizeboot+1:iboot*isizeboot,1)=ia(~isnan(ia));
    end
    xi=linspace(min(MatBoot),max(MatBoot),nbins);
    n_elements = histc(MatBoot,xi);
    c_elements = cumsum(n_elements);
    c_elements=c_elements./max(c_elements);
    %% the means for all windows
    pairs=a(~isnan(a));
    for ipairs=1:Nlink
        iTDM(ipairs,imethod)=c_elements(find(xi>=pairs(ipairs),1));
    end
    
    
end

%filec=load(datafile);

    Connectivity=zeros(Nnode);
MatC=Connectivity;
MatC=mln_removediag(MatC);

iTDM(:,Ninput+3)=MatC./max(MatC);
       
               
        