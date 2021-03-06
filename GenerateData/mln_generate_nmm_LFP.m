function mln_generate_nmm_LFP(is,nc,strfile,dirname,prename,cs,N,nl) 

% Huifang Wang For nmm Structure cluster Version
% related to the structures and pink noise
% Karl Friston
% $Id: spm_csd_demo.m 4402 2011-07-21 12:37:24Z karl $  
 dataname=[prename,'N',num2str(nc),'L',num2str(nl),'CS',num2str(100*cs),'S',num2str(is),'N',num2str(N)];
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
model        = 'LFP'; % or 'LFP' or 'CMC'
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
M.g   = 'spm_fx_lfp';
M.f   = f;
M.x   = x;
M.n   = numel(x(:));
M.m   = nc;
M.l   = nc;
 
 
 
% or generate data
%==========================================================================
 
% Integrate with pink noise process
%--------------------------------------------------------------------------
U.dt = 4/1000;
U.u  = randn(N,M.m)/100;
U.u  = sqrt(spm_Q(1/16,N))*U.u;
 
load(strfile,'PGS');
%===============================================================
P      = pE;
Astru  = PGS{is};
for  i = 1:size(Astru,2)
    P.A{1}(Astru(1,i),Astru(2,i)) = cs;
end
 
LFP     = spm_int_LFP(P,M,U);

LFP=LFP(100:N,:);
Data=full(LFP');
Params.fs=1000/8;
A1=P.A{1};
Connectivity=A1.*double(A1>=0);

 save(filename,'Data','Connectivity','Params');
 %msgbox('Data Generated','Success','custom',cdata);

 end

