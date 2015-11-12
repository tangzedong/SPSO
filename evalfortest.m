function [ objective ] = evalfortest( Task,rnvec,p_il,options )
%EVALFORTEST Summary of this function goes here
%   Detailed explanation goes here
d = Task.dims;
nvars = rnvec(1:d);
minrange = Task.Lb(1:d);
maxrange = Task.Ub(1:d);
y=maxrange-minrange;
vars = y.*nvars + minrange;
x=vars;
objective=Task.fnc(x);
end

