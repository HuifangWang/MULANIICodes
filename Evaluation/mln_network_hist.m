function fMat=mln_network_hist(Mat,bootn,thx,thy)
xvalues=linspace(0,1,bootn);
[Nnode,~,nboots]=size(Mat);
indx=find(xvalues>=thx,1);
fMat=NaN(Nnode);
Ny=thy*nboots;
for inode=1:Nnode
    for jnode=1:Nnode  
        pij=squeeze(Mat(inode,jnode,:));
        [yij,xij]=hist(pij,xvalues);
        if sum(yij(indx:end))>Ny
        fMat(inode,jnode)=1;
        else
            fMat(inode,jnode)=0;
        end
    end
end
