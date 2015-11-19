function [ obj ] = Perm( var,M,opt )
%AXIS Summary of this function goes here
%   Detailed explanation goes here
    dim = length(var);
%     opt=0*ones(1,dim);
%     opt(1:ceil(dim * 0.3)) = -10 * ones(1,ceil(dim * 0.3));
    var = (M*(var-opt)')';
    outtersum = 0;
    for k = 1:dim
        outtersum = outtersum + (sum( ([1:dim].^k + 0.5) .* ( (var(1:dim) ./ [1:dim]).^k - 1) ))^2;
    end
    obj = outtersum;
end