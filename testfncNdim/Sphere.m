function obj = Sphere(var,M,opt)
     D = length(var);
%     opt=0*ones(1,D);
    var = (var - opt);
    obj=var*var';
end