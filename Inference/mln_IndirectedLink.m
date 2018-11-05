%% Function to calculate the indirectedLink on A, 1 order is indirected link via (2^1) middle link, 2 order gets the link via at most 4 (2^2) links  
% Huifang Wang, Aug. 29 2013, Marseille
function C1=mln_IndirectedLink(A,Norder)

C1=firstorderindirectedlink(A);
if Norder>1
    for iorder=1:Norder
        B=C1+A;
        C1=firstorderindirectedlink(B);
        C1(C1.*A>0)=0;
    end
end


%% First order of the indirected link
function C1=firstorderindirectedlink(A)

Nchan=size(A,1);
C1=zeros(Nchan);
for irow=1:Nchan
idrow=find(A(irow,:));
idcol=find(A(:,irow));
C1(idcol,idrow)=1;
end
C1=C1-diag(diag(C1));
