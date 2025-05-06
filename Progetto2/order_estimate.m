function [p, C] = order_estimate(h, err)
% this function provides the asymptotic order of convergence and 
% the asymptotic constant of the sequence of errors contained in the vector err
%
% [p, C] = order_estimate(h, err)
%
p = log(err(2:end)./err(1:end-1))./log(h(2:end)./h(1:end-1));
C = err(2:end)./(h(2:end).^p);
%
return
