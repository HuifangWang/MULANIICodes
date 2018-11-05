function c=mln_Sim_Cos(M1,M2)
A=M1(:);
B=M2(:);
normA = sqrt(sum(A .^ 2));
normB = sqrt(sum(B .^ 2));

c = bsxfun(@rdivide, bsxfun(@rdivide, sum(A .* B), normA), normB);

if isnan(c)
    c=0;
end