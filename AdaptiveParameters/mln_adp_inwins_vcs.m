function outputname=mln_adp_inwins_vcs(dir,prenom,para,inwins,add,indtype,opt)
% this function is to output mulan for different numbers of windows

dir=[dir,'/ToutResults/'];

Ntheta=2;
if nargin<5
    add='';
end
[nop,npara]=size(opt.BestFPNch);
switch indtype
    case 'tak'
        if nop>=para.nboot2
        ind=randperm(nop,para.nboot2);
        
        else
            ind=repmat(1:nop,1,floor(para.nboot2/nop));
            Nmodeind=para.nboot2-length(ind);
            if Nmodeind>0
            ind=[ind,randperm(nop,Nmodeind)];
            end
        end
        
        Populations=opt.BestFPNch(ind,1:npara-Ntheta);
        para.theta=opt.BestFPNch(ind,npara-Ntheta+1:npara);
        
    case 'mean'
        pE=opt.mean;
        var=opt.var;
    tPopulations=mln_normrnd(pE,var,para.nboot2); % here nboots is the number of groups of parameters
    Populations=tPopulations(:,1:npara-Ntheta);
    para.theta=repmat(pE(npara-Ntheta+1:npara),[para.nboot2,1]);
    case 'hp'
        pE=opt.hp;
        var=opt.var;
    tPopulations=mln_normrnd(pE,var,para.nboot2); % here nboots is the number of groups of parameters
    Populations=tPopulations(:,1:npara-Ntheta);
    para.theta=repmat(pE(npara-Ntheta+1:npara),[para.nboot2,1]);
    
    case 'media'
        pE=opt.PE;
    var=opt.var;
    tPopulations=mln_normrnd(pE,var,para.nboot2); % here nboots is the number of groups of parameters
    Populations=tPopulations(:,1:npara-Ntheta);
    para.theta=repmat(pE(npara-Ntheta+1:npara),[para.nboot2,1]);
    
end
dataprenom=[dir,'Tout_',prenom];
outputname=['Adp_',prenom,add,'.mat'];
savestatisfile=[dir,outputname];
if ~exist(savestatisfile,'file')
     
    
    iMed=para.methodlog;
    Mati=load(dataprenom,iMed{:});
    %% given interested Nwin
    
    [Nchan,~, Nwins]=size(Mati.(iMed{1}));
    
    if nargin<4 || inwins==-1
        Setlist{1}=1:Nwins;
    else
        Setlist=mln_set_winlist(Nwins,inwins);
    end
    
    %%
    mlnMat=NaN(Nchan,Nchan,para.nboot2,length(Setlist));
   
    para.Setlist=Setlist;
    % read all used methods
    for ipop=1:para.nboot2
        para.pc = Populations (ipop,:);
        for iwins=1:length(Setlist)
            methodlog=para.methodlog;
            MatA=load(dataprenom,methodlog{:});
            InputMM=mln_GenerateTDM_win_dy(MatA,methodlog,para.nboot,1,Setlist{iwins});
            
              
                        
            %Ninput=length(methodlog);
            pause(0.01)
            indMM=InputMM.methods;
            %ind_methods=[1 3 5 6]';%
            ind_methods=mln_chs_ind(indMM,methodlog);
            iInput=InputMM.values(:,ind_methods);
               
            [iMat,~] = mln_algorithm_para(iInput,para,para.kforMF,para.methodlog);
            mlnMat(:,:,ipop,iwins)=iMat;

                        
         end
               
                
     end
            
        

    para.Population=Populations;
    save(savestatisfile,'mlnMat','para');
end


function Input=mln_GenerateTDM_win_dy(MatA,methods_used,nboot,nboot2,winsV)

Nchan=size(MatA.(methods_used{1}),1);
Nlink=Nchan*(Nchan-1);
Nmethod=length(methods_used);
Input.values=NaN(Nlink,Nmethod,nboot2);
Input.methods=methods_used;

nbins=100;
Norder=2;

for imethod=1:Nmethod
    Mat1=MatA.(methods_used{imethod});
    Mat1=Mat1(:,:,winsV);
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

 
        
function Setlist=mln_set_winlist(Nwins,inwins)

Nmln=Nwins/inwins;
for imln=1:Nmln
    Setlist{imln}=(imln-1)*inwins+1:imln*inwins;
end


function Population=mln_normrnd(pE,var,nboot2)

Population=zeros(nboot2,length(pE));
Population(1,:) = pE;
for iPop = 2:nboot2
    Population(iPop,:) = normrnd(pE,var);
end
