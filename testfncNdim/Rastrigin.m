function obj = Rastrigin(var,M,opt)
    dim = length(var);
%     opt=0*ones(1,dim);
%     opt=49*ones(1,dim);
%     opt(1:ceil(dim * 0.3)) = 20 * ones(1,ceil(dim * 0.3));
    var = (M*(var-opt)')';
    obj = 10*dim;
    for i=1:dim
        obj=obj+(var(i)^2 - 10*(cos(2*pi*var(i))));
    end
end