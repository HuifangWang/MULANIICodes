function mln_Gen_NMM_fs_vs_ns_few(dirname,prename,nc,is,npts,mincs,maxcs,nl,fs)
% generate the data for given fs with varying cs

% dirname: the folder to store all files
% prenom: spectial name for new datasets
% nc: the channels of data
% is: which structure stored in the structure file
% npts: length of data
% cs: the connection strength
% nl: how many strong links

% -----------------------------------------------------
% Examples:mln_Gen_NMM_fs_vs Examples nmm 20 1 810 0.1 1 10 250
% Note that: there should be a structureN*L*S10.mat file in the dirname
% folder
% structureN*L*S10.mat can be generated by calling function
% './FindStructures/mln_gen_structure_Nodes.m'
% -------------------------------------------------------------------------

% Huifang Wang, Marseille, Dec 3, 2014
% -------------------------------------------------------------------------

datadir=[dirname,'/data'];
if ~exist(datadir,'dir')
    try mkdir(datadir);
    end
end


is=str2double(is);
nc=str2double(nc);
nl=str2double(nl);

npts=str2double(npts);
mincs=str2double(mincs);
maxcs=str2double(maxcs);
fs=str2double(fs);   


mln_generate_nmm_vcs_fs_few(is,nc,dirname,prename,mincs,maxcs,npts,nl,fs); 