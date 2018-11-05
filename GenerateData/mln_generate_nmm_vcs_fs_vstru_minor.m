% Huifang Wang For nmm Structure cluster Version
% related to the structures and pink noise
% Karl Friston
% $Id: spm_csd_demo.m 4402 2011-07-21 12:37:24Z karl $
% Feb 13, Huifang Wang dynamic
%-------------------------------------------------------------
%Nt: the time periods for Nstru structures
 function dataname=mln_generate_nmm_vcs_fs_vstru_minor(is,nc,strfile,dirname,prename,mincs,maxcs,Nt,nl,fs,Nstru,chanp)

 ncl=floor((nl*nc/5)*chanp);
 dataname=[prename,'N',num2str(nc),'L',num2str(nl),'CS',num2str(100*mincs),num2str(100*maxcs),'S',num2str(is),'Nt',num2str(Nt),'Nstru',num2str(Nstru),'MiCNL',num2str(ncl)];
 filename=['./',dirname,'/data/',dataname,'.mat'];
 %[cdata] = imread('MULANLOGO.png'); 
 if exist(filename,'file')
     disp('Data already there')
     %msgbox('Data already there','Success');
 else
     if ~exist(dirname,'dir')
         mkdir(dirname);
     end
     datadir=[dirname,'/data'];
     if ~exist(datadir,'dir')
         mkdir(datadir);
     end

% number of sources and LFP channels (usually the same)
%--------------------------------------------------------------------------
model        = 'ERP'; % or 'LFP' or 'CMC'
dipfit.type  = 'LFP';
dipfit.model = model;
dipfit.Ns    = nc;
dipfit.Nc    = nc;
 
 
% specify network (connections)
%--------------------------------------------------------------------------
A{1}  = zeros(nc);                            % a forward connection
A{2}  = zeros(nc);                            % a backward connection
A{3}  = sparse(nc,nc);                         % lateral connections
B     = {};                                  % trial-specific modulation
C     = speye(nc,nc);                          % sources receiving innovations
 
% get priors
%--------------------------------------------------------------------------
pE    = spm_dcm_neural_priors(A,B,C,dipfit.model);  % neuronal priors
pE    = spm_L_priors(dipfit,pE);          % spatial  priors
[x,f] = spm_dcm_x_neural(pE,dipfit.model);
 
% create LFP model
%--------------------------------------------------------------------------
M.dipfit = dipfit;
M.g   = 'spm_gx_erp';
M.f   = f;
M.x   = x;
M.n   = numel(x(:));
M.m   = nc;
M.l   = nc;
 
preN=ceil(10*fs);%second 
 
% or generate data
%==========================================================================
 
% Integrate with pink noise process
%--------------------------------------------------------------------------
U.dt = 1/fs;


U.u  = randn(preN,M.m)/100;
U.u  = sqrt(spm_Q(1/16,preN))*U.u;
 
load(strfile,'PGS');
%===============================================================

P      = pE;
A10=P.A{1};
%Vstru=[1,1,1,1];
Astru  = PGS{is};
minorPGS=mln_upd_minorPGS(Astru,ncl,Nstru,nc);% ncl is the number of links which will be changed
P.A{1}=mln_updateP(A10,Astru,mincs,maxcs);
%%% generate the first stable signals
[~, x0] = spm_int_L(P,M,U);
%M.x=spm_unvec(x0,M.x);
iNtp=ceil(Nt*fs);
Data=NaN(nc,Nstru*iNtp);
Connects=NaN(nc,nc,Nstru);
for istr=1:Nstru
U.u  = randn(iNtp,M.m)/100;
U.u  = sqrt(spm_Q(1/16,iNtp))*U.u;    
P.A{1}=mln_updateP(A10,minorPGS{istr},mincs,maxcs);
[LFP,x0] = spm_int_Lx0(x0,P,M,U);
%M.x=spm_unvec(x0,M.x);
Data(:,(istr-1)*iNtp+1:istr*iNtp)=LFP';
A1=P.A{1};
Connects(:,:,istr)=A1.*double(A1>=0);
end
Connectivity.strus=Connects;
Params.fs=fs;
Connectivity.switchTimes=(0:Nstru-1)*Nt;

save(filename,'Data','Connectivity','Params');
%msgbox('Data Generated','Success','custom',cdata);

end
 

function CA=mln_updateP(CA,Astru,mincs,maxcs)
Ncontcs=10;
for  i = 1:size(Astru,2)
    ics=randi([floor(mincs*Ncontcs) floor(maxcs*Ncontcs)])/Ncontcs;
    CA(Astru(1,i),Astru(2,i)) = ics;
end

function minorPGS=mln_upd_minorPGS(Astru,ncl,Nstru,nc)
Nlinks=size(Astru,2);
minorPGS{1}=Astru;
for istru=2:Nstru
    chanlink=randperm(Nlinks,ncl);
    oldStru=minorPGS{istru-1};
    newStru=oldStru;
    for incl=1:ncl
        newlink=mln_newlinkp1(oldStru(2,chanlink(incl)),nc);
        if newlink==oldStru(1,chanlink(incl))
            newlink=mln_newlinkp1(newlink,nc);
        end
        newStru(2,chanlink(incl))=newlink;
    end
    minorPGS{istru}=newStru;
end

function newlink=mln_newlinkp1(newlink,nc)
newlink=newlink+1;
if newlink>nc
    newlink=1;
end
    



