function [x, niter, err] = fixed_point (x0, T, toll, itmax)
% Author: R. Sacco
%
% This function implements the fixed point iteration method 
% to numerically compute the zero of g=g(x)
% 
% [x, niter, err] = fixed_point (x0, T, toll, itmax)
%
x    = x0;
xold = x;
%
k    = 0;
e    = toll + 1;
err  = [];
%
while ((k < itmax) & (e > toll))
      %
      xnew = T(xold);
      %
      if (isnan(xnew) | isinf(xnew))
         %
         e = 0;
         %
      else
         %
         incr = xnew - xold;
         e    = abs(incr);
         x    = [x; xnew];
         err  = [err; incr];
         xold = xnew;
         k    = k+1;
         %
      end
end
%
niter = k;
%
return
