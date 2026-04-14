function [p1, p2, a, b] = fastContactsolverTesting(param, t1, t2, z0)
% Wrapper that prefers the compiled MEX if available (same API).
% Default behavior: uses the same initial guess for every sample (no warm-start).
%
% Need to run these lines of code in command prmpt first:
% clear fastContactSolver
% delete(gcp('nocreate'));
% mex -O fastContactSolver_mex.c
% which -all fastContactSolver_mex
%
% To revert to your pure MATLAB implementation for testing, rename this file
% or temporarily shadow fastContactsolver_mex.

if nargin < 4, z0 = []; end
if exist('fastContactSolver_mex','file') == 3
    if isempty(z0), [p1,p2,a,b] = fastContactSolver_mex(param, t1, t2);
    else,           [p1,p2,a,b] = fastContactSolver_mex(param, t1, t2, z0);
    end
else
    % Fallback: call your original MATLAB code (put it in a private function
    % or separate file named fastContactsolver_matlab.m).
    [p1,p2,a,b] = fastContactSolver(param, t1, t2, z0);
end
end