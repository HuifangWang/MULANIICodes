 % Huifang Wang For nmm Structure cluster Version
% related to the structures and pink noise
% Karl Friston
% $Id: spm_csd_demo.m 4402 2011-07-21 12:37:24Z karl $
 function dataname=mln_generate_erp_fs_hc(is,nc,strfile,dirname,prename,N,fs) 
  
 dataname=[prename,'N',num2str(nc),'S',num2str(is),'N',num2str(N)];
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

%===============================================================
P      = pE;


load(strfile,'sMat');
isMat=sMat(:,:,is);
%===============================================================
%minMat=min(min(P.A{1}));
%P.A{1} = isMat;

%P.A{1}(P.A{1}<0.5) = minMat;

for inch= 1:nc
    for jinch = 1:nc
        ics =isMat(inch,jinch);
        if ics>=0.06
            P.A{1}(inch,jinch)=ics;
        end
        
    end
end
U.dt = 1/fs;

for icount=1:10

% Integrate with pink noise process
%--------------------------------------------------------------------------

U.u  = randn(N,M.m)/100;
U.u  = sqrt(spm_Q(1/16,N))*U.u;
 
LFP     = spm_int_L(P,M,U);
LFP=LFP(100:N,:);
Data=full(LFP');
if max(max(Data))-min(min(Data))<0.0035
    break;
end
end
if max(max(Data))-min(min(Data))>0.0035
    return;
end
Params.fs=fs;
A1=P.A{1};
Connectivity=A1.*double(A1>=0);

 save(filename,'Data','Connectivity','Params');
 %msgbox('Data Generated','Success','custom',cdata);

 end
