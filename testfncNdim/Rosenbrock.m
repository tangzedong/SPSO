function [ obj ] = Rosenbrock( var,M,opt )
%AXIS Summary of this function goes here
%   Detailed explanation goes here
    dim = length(var);
%     opt=0*ones(1,dim);
%     opt(1:ceil(dim * 0.3)) = -10 * ones(1,ceil(dim * 0.3));
    var = (M*(var-opt)')';
    obj = sum((var(2:end) - var(1:end-1).^2).^2 + (1-var(1:end-1)).^2);
end