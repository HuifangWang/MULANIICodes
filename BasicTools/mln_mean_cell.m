
function mM=mln_mean_cell(cM)
dim = ndims(cM{1});          
M = cat(dim+1,cM{:});        %# Convert to a (dim+1)-dimensional matrix
mM = mean(M,dim+1);