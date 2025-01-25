function CoverTree = covertree_build( X, opts )

%
% function [idxs,dists] = covertree_build( cX, opts )
%
% Wrapper for cover tree covertree construction function
%
% IN:
%   cX              : D by N matrix of N points in R^D
%   [opts]          : structure of options.
%                       [NTHREADS]      : number of threads for search. Default: number of cores.
%
% OUT:
%   CoverTree       : covertree data structure, with the following fields:
%                       theta   : dilation factor (0.5 is the default)
%                       radii   : row vector of radii at different levels
%                       levels  : n by 5 matrix. the i-th row contains [level,parent,#children,index first children,children idxs]
%                                 for the i-th point.
%                       outparams : 9 vector of parameters
%                       nclasstogetdist : number of distance computations
%
%
% (c) Copyright Duke University, 2014
% Mauro Maggioni
% mauro@math.duke.edu
%
persistent NTHREADS

CoverTree = [];

MIN_N_FOR_MULTITHREADED  = 4096;

if nargin<2,                    opts = struct();                            end
if size(X,2)>MIN_N_FOR_MULTITHREADED
    if isfield(opts,'NTHREADS'),
        NTHREADS = int32(opts.NTHREADS);
    else
        NTHREADS = int32(feature('numcores'));
    end
else
    NTHREADS = int32(0);
end
if ~isfield(opts,'BLOCKSIZE'),  opts.BLOCKSIZE = int32(2048);                end

opts.NTHREADS       = NTHREADS;
opts.theta          = single(opts.theta);
opts.distancefcn    = int32(0);
opts.classname      = int32(0);
opts.numlevels      = int32(opts.numlevels);

if ~isempty(X)
    if isa(X,'single')
        CoverTree = covertree( opts, X );
    elseif isa(X,'double')
        opts.theta    = double(opts.theta);
        CoverTree = covertreeD( opts, X );
    else
        warning('\n covertree_build: invalid class %s for X',class(X));
    end
end

return