%% First cirterion to stop the MIS itinerary, the same structure; 
%% Huifang Wang 24/sep./2013, Marseille
function isstop=mln_Miscriteria(M1,M2)
Nnode=size(M1,1);
B1=M1~=0;
B2=M2~=0;
A=B1==B2;
if length(find(A==1))==Nnode*Nnode
    isstop=1;
else
    isstop=0;
end


