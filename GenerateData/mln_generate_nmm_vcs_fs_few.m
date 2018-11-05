 % Huifang Wang For nmm Structure cluster Version
% related to the structures and pink noise
% Karl Friston
% $Id: spm_csd_demo.m 4402 2011-07-21 12:37:24Z karl $
 function dataname=mln_generate_nmm_vcs_fs_few(is,nc,dirname,prename,mincs,maxcs,N,nl,fs) 
  
 dataname=[prename,'N',num2str(nc),'L',num2str(nl),'CS',num2str(100*mincs),num2str(100*maxcs),'S',num2str(is),'N',num2str(N)];
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
 
 
 
% or generate data
%==========================================================================
Ntimes=10; 
% Integrate with pink noise process
%--------------------------------------------------------------------------
U.dt = 1/fs;
U.u  = randn(N,M.m)/100;
U.u  = sqrt(spm_Q(1/16,N))*U.u;
 
%===============================================================
P      = pE;
%%% nl determines how many links from 0.8 to maxcs
%%% nlw determines how many links from 0.5 to 0.8
%%% nln determines how many links from 0.1 to 0.5
%nl=20;
nlw=nl;
nln=nl;
nll=[nl,nlw,nln];
mincs=[0.8,0.5,0.1];
maxcs=[1,0.7,0.4];
for ilevel=1:length(nll)
    for inl = 1:nll(ilevel)
        ics=randi([floor(mincs(ilevel)*Ntimes) floor(maxcs(ilevel)*Ntimes)])/Ntimes;
        
        linkft=randi(nc,[2,1]);
        while linkft(1)==linkft(2)
          linkft=randi(nc,[2,1]);
        end
        P.A{1}(linkft(1),linkft(2)) = ics;
        
    end
end
 
LFP     = spm_int_L(P,M,U);
LFP=LFP(100:N,:);
Data=full(LFP');
Params.fs=fs;
A1=P.A{1};
Connectivity=A1.*double(A1>=0);

 save(filename,'Data','Connectivity','Params');
 %msgbox('Data Generated','Success','custom',cdata);

 end

