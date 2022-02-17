function varargout = calcmld(depth, sigtemp, cutoff, refdepth)
%CALCMLD Calculate mixed layer depth
%
% mld = calcmld(depth, temp)
% mld = calcmld(depth, sig)
% mld = calcmld(..., cutoff, refdepth)
% [mld, mldmask] = calcmld(...)
%
% Input variables:
%
%   depth:      vector of length nd corresponding to depths of sig/temp
%               grid 
%
%   sig:        nd x m array of sigma-t values
%
%   temp:       nd x m array of temperature values
%
%   cutoff:     Cutoff gradient used to define mixed layer.  Default is
%               0.125 for density input, and 0.5 for temperature input
%               (Levitus 1982)
%
%   refdepth:   Reference depth used to compare gradient to.  Default is
%               surface (i.e. shallowest depth available)
%
% Output variables"
%
%   mld:        1 x m vector of mixed layer depths
%
%   mldmask:    nd x m logical array indicating whether each input value is
%               in the mixed layer (true) or not (false)

% Copyright 2009 Kelly Kearney

%------------------------
% Parse input
%------------------------

if ~isvector(depth)
    error('Depth must be vector');
end

depth = depth(:);
nd = length(depth);

sz = size(sigtemp);
if ~any(sz == nd)
    error('Density or temperature array does not match length of depth vector');
end

if sz(1) ~= nd && sz(2) == nd
    sigtemp = sigtemp';
end
nt = size(sigtemp,2);

if nargin < 3 || isempty(cutoff)
    if all(sigtemp > 1000)  % Density
        cutoff = 0.125;
    else
        cutoff = 0.5;       % Temperature
    end
end

if nargin < 4 || isempty(refdepth)
    refidx = 1;
else
    [blah, refidx] = min(abs(depth - refdepth));
end


%------------------------
% Calculate mixed layer
% depth
%------------------------

mld = zeros(nt,1);

sigdiff = abs(bsxfun(@minus, sigtemp, sigtemp(1,:)));

for it = 1:nt
    isbad = isnan(sigdiff(:,it));
    x = sigdiff(~isbad, it);
    y = depth(~isbad);
    
    idx = find(x > cutoff, 1, 'first');
    if isempty(idx)
        mld(it) = NaN;
    else
        mld(it) = interp1q(x(idx-1:idx), y(idx-1:idx), cutoff);
    end
    
end

if nargout == 1
    varargout{1} = mld;
elseif nargout == 2
    varargout{1} = mld;
    varargout{2} = sigdiff < cutoff;
end
    