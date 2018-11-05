function P=mln_IcovPartial(M,lambda)

%% This function transfers a Matrix M to Partial ones P
%lambda=5,100;
nchannel=size(M,1);
oc=M-diag(diag(M))+eye(nchannel);
IM=-L1precisionBCD(oc/mean(diag(oc)),lambda/1000);
P=(IM ./ repmat(sqrt(abs(diag(IM))),1,nchannel)) ./ repmat(sqrt(abs(diag(IM)))',nchannel,1)+eye(nchannel);


