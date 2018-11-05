function mf_param = mln_genparam(data,mf_n,mf_type,k)
%GENPARAM Generate initial membership function parameters for ANFIS learning.
%   GENPARAM(DATA,MF_N,MF_TYPE) generates initial input MF parameters
%   from a M-by-N training data matrix DATA, where M is the number of
%   training data pairs and N is the number of inputs plus one.
%   MF_N and MF_TYPE are optional arguments that specify the MF
%   number and MF type of each input, respectively. MF_N should be a
%   vector of length N. If MF_N is a number, then it applies to all inputs.
%   Similarly, MF_TYPE should be a string matrix of N rows. If MF_TYPE is a
%   single string, it applies to all inputs. Default values for MF_N and
%   MF_TYPE are 2 and 'gbellmf', respectively.
%
%   The centers of the generated MFs are always equally spaced along the
%   domain of an input variable, where the domain is determined as the
%   interval between the min. and max. of the corresponding column in DATA.
%
%   Restrictions: (1) 'sigmf', 'smf', and 'zmf' MF types are not supported
%   since they are either open-left or open-right. (2) The same MF type
%   is assigned to MFs of the same input variable.
%
%   For example:
%
%       NumData = 1000;
%       data = [rand(NumData,1) 10*rand(NumData,1)-5 rand(NumData,1)];
%       NumMf = [3 7];
%       MfType = str2mat('trapmf','gbellmf');
%       MfParams = genparam(data,NumMf,MfType);
%       set(gcf,'Name','genparam','NumberTitle','off');
%       NumInput = size(data, 2) - 1;
%       range = [min(data)' max(data)'];
%       FirstIndex = [0 cumsum(NumMf)];
%       for i = 1:NumInput;
%           subplot(NumInput, 1, i);
%           x = linspace(range(i, 1), range(i, 2), 100);
%           index = FirstIndex(i)+1:FirstIndex(i)+NumMf(i);
%           mf = evalmmf(x, MfParams(index, :), MfType(i,:));
%           plot(x, mf');
%           xlabel(['input ' num2str(i) ' (' MfType(i, :) ')']);
%       end
%
%   See also GENFIS1, ANFIS.

%       Roger Jang, 8-7-94
% Huifang Wang, Sep. 4, 2013
%   Copyright 1994-2006 The MathWorks, Inc.
%   $Revision: 1.14.2.1 $  $Date: 2006/09/30 00:18:46 $

% Change this to have different default values
default_mf_n = 2;
default_mf_type = 'gbellmf';

if nargin <= 2,
    mf_type = default_mf_type;
end
if nargin <= 1,
    mf_n = default_mf_n;
end

if nargin <= 3,
    k=5;
end
in_n=length(data);

% generate mf_n and mf_type of proper sizes
if length(mf_n) == 1,
    mf_n = mf_n(ones(in_n, 1), :);
end
if size(mf_type, 1) == 1,   % single mf_type for all MFs
    mf_type = mf_type(ones(in_n, 1), :);
end


mf_param = [];
for i = 1:in_n;
    type = deblank(mf_type(i, :));
    if strcmp(type, 'gbellmf'),
    
        a = [data(i),1-data(i)]';
        b = k*2.*a;
        c =[0,1]';
        
        mf_param = [mf_param; 
            a b c];
    elseif strcmp(type, 'gaussmf'),
        a = (range(i,2) - range(i, 1))/2/(mf_n(i) - 1);
        a = a/sqrt(2*log(2));
        c = linspace(range(i, 1), range(i, 2), mf_n(i))';
        d = zeros(mf_n(i), 1);
        mf_param = [mf_param; a(ones(mf_n(i), 1)) c d d];
    elseif strcmp(type, 'gauss2mf'),
        a = (range(i,2) - range(i, 1))/2/(mf_n(i) - 1);
        c = linspace(range(i, 1), range(i, 2), mf_n(i))';
        c1 = c - 0.6*a;
        c2 = c + 0.6*a;
        sigma = 0.4*a/sqrt(2*log(2))*ones(mf_n(i), 1); 
        mf_param = [mf_param; sigma c1 sigma c2];
    else
        fprintf('mf_type = %s\n', type);
        error('Unsupported MF type!');
    end
end