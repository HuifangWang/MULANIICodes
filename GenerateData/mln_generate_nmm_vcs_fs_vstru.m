% Huifang Wang For nmm Structure cluster Version
%dd
% related to the structures and pink noise
% Karl Friston
% $Id: spm_csd_demo.m 4402 2011-07-21 12:37:24Z karl $
% Feb 13, Huifang Wang dynamic
%-------------------------------------------------------------
%Nt: the time periods for Nstru structures
 function dataname=mln_generate_nmm_vcs_fs_vstru(is,nc,strfile,dirname,prename,mincs,maxcs,Nt,nl,fs,Nstru)

  
 dataname=[prename,'N',num2str(nc),'L',num2str(nl),'CS',num2str(100*mincs),num2str(100*maxcs),'S',num2str(is),'Nt',num2str(Nt),'Nstru',num2str(Nstru)];
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
NPGS=length(PGS);
randint=randperm(NPGS);
Vstru=randint(1:Nstru);
%Vstru=[1,1,1,1];
Astru  = PGS{Vstru(1)};
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
P.A{1}=mln_updateP(A10,PGS{Vstru(istr)},mincs,maxcs);
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


