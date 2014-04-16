function [y] = digamma(x)

% function [y] = digamma(x)
%
% The digamma function (derivative of the gamma function)
%
% Implementation via MEX routine
%
% Solaris: digamma.mexsol
%
% Compiled via make function (based on mex digamma.c)


y = psi(x);