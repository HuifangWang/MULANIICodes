function AUCthmmerge=mln_merge_ind(M,indM,merge_ind,rest_ind)

%Nrd=length(rest_ind);
[indmd,~]=mln_chs_ind(indM,merge_ind);
[indrd,~]=mln_chs_ind(indM,rest_ind);
%indMnew=[merge_ind,rest_ind];
indswap=[indmd',indrd'];
M=permute(M,indswap);
indmd=1:length(indmd);
indrd=length(indmd)+1:length(indmd)+length(indrd);



sizeM=size(M);
AUCthmmerge=cell(sizeM(indrd));




   alld=numel(AUCthmmerge);

  Ncolon=length(indmd);
  indcolon=repmat(':,',1,Ncolon-1); 
   indcolon=[indcolon,':,'];
         

   
   for id=1:alld
       AUCthmmerge{id}=eval(['M(',indcolon,'id)']);
   end
   
   