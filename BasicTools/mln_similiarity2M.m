function S2M=mln_similiarity2M(A,B)

A1=~~A;
B1=~~B;
[Nchan,~]=size(A);
Nlink=Nchan*(Nchan-1);
A1=A1-diag(diag(A1));
B1=B1-diag(diag(B1));
AB=A1==B1;

S2M=(length(find(AB==1))-Nchan)/Nlink;

