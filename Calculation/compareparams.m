function issame=compareparams(oldparams,newparams)
%% compare two sets of parameters
%fieldold=fieldnames(oldparams);
fieldnew=fieldnames(newparams);
Nfield=length(fieldnew);
issame=0;
for i=1:Nfield
    ifield=fieldnew{i};
    %Nifield=length(ifield);
    %find(strncmpi(fileprevar,ifield,Nifield)==1);
    if oldparams.(ifield)==newparams.(ifield);
        issame=1;
    else
        issame=0;
        return;
    end
end
