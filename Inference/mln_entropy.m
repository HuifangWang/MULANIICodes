function  H = mln_entropy(X,nb)
%xr=linspace(min(X),max(X),nb);
xr=linspace(0,1,nb);
p = hist(X,xr);
dp=xr(2)-xr(1);

nz=p>0;
freq = p(nz)/sum(p);
H = -sum(freq.*dp.*log2(freq));

% freq =p(nz)/sum(p);
% H = -sum(freq.*dp.*log(freq));

%H = -sum(freq.*log2(freq));

