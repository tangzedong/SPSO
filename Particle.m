classdef Particle    
    properties
        uvec; % (genotype)--> decode to find design variables --> (phenotype) 
        velocity; %individial velocity
        lbvec; %individial best vector
        
        %individial best feature
        lfactorial_costs;
        lfactorial_ranks;
        lscalar_fitness;
        lskill_factor;
        
        %current particle feature
        factorial_costs;
        factorial_ranks;
        scalar_fitness;
        skill_factor;
    end    
    methods        
        function object = initialize(object,D)            
            object.uvec = rand(1,D);
            object.velocity = 0.1 * rand(1, D);
            object.lbvec = object.uvec;
        end
        
        function [object,calls] = evaluate(object,Tasks,p_il,no_of_tasks,options)     
            if object.skill_factor == 0
                calls=0;
                for i = 1:no_of_tasks
                    [object.factorial_costs(i),xxx,funcCount]=fnceval(Tasks(i),object.uvec,p_il,options);
                    calls = calls + funcCount;
                end
            else
                object.factorial_costs(1:no_of_tasks)=inf;
                for i = 1:no_of_tasks
                    if object.skill_factor == i
                        [object.factorial_costs(object.skill_factor),object.uvec,funcCount]=fnceval(Tasks(object.skill_factor),object.uvec,p_il,options);
                        calls = funcCount;
                        break;
                    end
                end
            end
        end
        
        function [object,calls] = evaluate_SOO(object,Task,p_il,options)   
            [object.factorial_costs,object.uvec,funcCount]=fnceval(Task,object.uvec,p_il,options);
            calls = funcCount;
        end
        
        function object=forward(object,p1,p2,cf)
            object.uvec=0.5*((1+cf).*p1.uvec + (1-cf).*p2.uvec);
            object.uvec(object.uvec>1)=1;
            object.uvec(object.uvec<0)=0;
        end  
    end
end