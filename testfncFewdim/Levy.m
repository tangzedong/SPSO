function [ obj ] = Levy( var,M,opt )
%AXIS Summary of this function goes here
%   Detailed explanation goes here
    dim = length(var);
%     opt=0*ones(1,dim);
%     opt(1:ceil(dim * 0.3)) = -10 * ones(1,ceil(dim * 0.3));
    var = (M*(var-opt)')';
    obj = sum([1:dim].*var(:));
end