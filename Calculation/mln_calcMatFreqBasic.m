function mln_calcMatFreqBasic(Resultfile,VMethlog,lfp,params)
%% to calculate the connectivity matrix for FreqBasic methods in Methlog

% Huifang Wang, June 9, 2012, Inserm U1106, Marseille
Methlog=char(VMethlog{1});
NMlog=length(Methlog);
if exist(Resultfile,'file')

        return;

end
%% if there is not the result then calculate

[Nchannel, Ntime]=size(lfp);
freqs=params.minfreq:params.stepfreq:params.maxfreq;
Nfreq=length(freqs);

if ~params.overlap==0
overlap_p=params.wins*params.overlap;
Nwindows=floor((Ntime-overlap_p)/(params.wins-overlap_p));
else
    overlap_p=0;
    Nwindows=floor(Ntime/params.wins);
end
% initialization
Dmatf=[Nchannel,Nchannel,Nfreq,Nwindows];
NMeths=length(VMethlog);
for i=1:NMeths
    iMeth=char(VMethlog{i});
    if ~isempty(iMeth)
Mat.(iMeth)=zeros(Dmatf); 
    end
end
VMethodlog=fieldnames(Mat);
Nmethod=length(VMethodlog);
for i=1:Nwindows
    i_lfp=lfp(:,floor((i-1)*(params.wins-overlap_p)+1):floor(i*params.wins-(i-1)*overlap_p));
    iMat=mln_icalcMatFreqBasic(i_lfp,freqs,params.fs);
    for j=1:Nmethod
            Mat.(char(VMethodlog(j)))(:,:,:,i)=iMat.(char(VMethodlog(j)));
    end
%     Mat.psd.w(:,floor((i-1)*params.wins+1):floor(i*params.wins),:)=iMat.cfs.w;
%     Mat.psd.f(:,i,:)=iMat.cfs.f;
end
updateResult(Resultfile,Mat,params);