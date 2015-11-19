function [data_PSO, ppt] = PSO(Tasks, option)
    clc    
    tic   
%%%初始化参数
    data_PSO.success = 1;
    c1 = option.c1;
    c2 = option.c2;
    momentum = option.momentum;
    pop = option.pop;
    maxgen = option.maxgen;
    p_il = option.p_il;
    rho = option.alpha;
    ntotalcall = 1;
    mu = 0.01;
    if mod(pop,2) ~= 0  %保证种群数是偶数
        pop = pop + 1;
    end
    ntasks = 1;
    dim = Tasks(1).dims;  %取任务中最大的维数作为个体向量的维数
    local_search_opt = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton');   %个体实施局部优化
    
    fnceval_calls = 0;      %计算目标函数执行次数
    calls_per_individual = zeros(1,pop);
    EvBestFitness = zeros(ntasks,maxgen);  %记录最好适应度
    TotalEvaluations = zeros(1,maxgen);
    bestobj = inf;      %记录找到的最优目标函数值
    funvalue = -1*ones(ntasks,pop);
%%%初始化种群
    for i = 1 : pop
        swarm(i) = Particle();         
        swarm(i) = initialize(swarm(i),dim);
    end
    for i = 1 : pop
        [swarm(i),calls_per_individual(i)] = evaluate_SOO(swarm(i),Tasks,p_il,local_search_opt);
        funvalue(ntotalcall) = swarm(i).factorial_costs;
        ntotalcall = ntotalcall + 1;
    end
    
    fnceval_calls = fnceval_calls + sum(calls_per_individual);        %累加每个个体的目标函数计算次数
    TotalEvaluations(1) = fnceval_calls; %TotalEvaluations记录每次迭代完成时目标函数的执行次数
    
%%%初始化种群的factorial_cost, factorial_ranks
    factorial_cost = zeros(1,pop);
    for j = 1:pop
        factorial_cost(j) = swarm(j).factorial_costs;
    end
    [xxx,y]=sort(factorial_cost);
    swarm = swarm(y);
    for j=1:pop
        swarm(j).factorial_ranks(i)=j;
    end
    bestobj = swarm(1).factorial_costs;   %bestobj(i)记录目标函数(i)的最好值
    oldbestobj = bestobj;
    EvBestFitness = bestobj;              %EvBestFitness(i)记录每一代的(i)最优值
    bestInd_data = swarm(1);              %bestInd_data(i)记录最优个体
    bestvec = swarm(1).uvec;
    %         bestfitness(i) = swarm(1).scalar_fitness;

    
    generation=0;
    eps = 0.1;
    while generation <= maxgen %&& bestobj(1) > eps
        generation = generation + 1;
        Velmax = 0.3 * exp(-20*(generation-1)/maxgen) + 0.006;
        disp(Velmax);
        nexchange = 0;
        for i = 1:pop
            swarm(i).success_exchange = 0;
        end
        for i = 1:pop
            factcost(i, :) = swarm(i).factorial_costs;
            velevery(i, :) = swarm(i).velocity;
            uvec(i, :) = swarm(i).uvec;
        end
%%to do 
%%pso的个体速度和位置的更新，选择最优个体，记录局部最优位置

        u = rand(1,dim);
        cf = zeros(1,dim);
        cf(u<=0.5)=(2*u(u<=0.5)).^(1/(mu+1));
        cf(u>0.5)=(2*(1-u(u>0.5))).^(-1/(mu+1));
        for i = 1 : pop
            Vel = swarm(i).velocity;
            Vel = momentum * Vel ...
                + c2 .* rand(1,dim) .* bsxfun(@minus, bestvec, swarm(i).uvec) ...
                + c1 .* rand(1,dim) .* bsxfun(@minus, swarm(i).lbvec, swarm(i).uvec);

            swarm(i).velocity = Vel;
            
            swarm(i).velocity(Vel < -Velmax) = -Velmax;
            swarm(i).velocity(Vel > Velmax) = Velmax;

            PS = swarm(i).uvec;

            PS = rho * PS + Vel;
            swarm(i).uvec = PS;
            swarm(i).uvec(PS < 0) = 0;
            swarm(i).uvec(PS > 1) = 1;
        end
        
        %%评估更新后的个体
        parfor a = 1 : pop
            [swarm(a),calls_per_individual(a)] = evaluate_SOO(swarm(a),Tasks,p_il,local_search_opt);
        end
        for i = 1:pop
            funvalue(ntotalcall) = swarm(i).factorial_costs;
            ntotalcall = ntotalcall + 1;
        end
        fnceval_calls = fnceval_calls + sum(calls_per_individual);
        TotalEvaluations(generation) = fnceval_calls;
        
        factorial_cost = zeros(1,pop);
        for j = 1:pop
            factorial_cost(j) = swarm(j).factorial_costs;
        end
        [~,ind] = sort(factorial_cost);
        swarm = swarm(ind);
        
        if swarm(1).factorial_costs <= bestobj
            bestvec(:) = swarm(1).uvec;
            bestobj = swarm(1).factorial_costs;
            bestInd_data(i) = swarm(1);
            data_PSO.bestindx = ind(1);
            data_PSO.nfc = ind(1) + generation * pop;
        end
        EvBestFitness(generation) = bestobj;
     
        for i=1:pop
            %记录个体历史最优值
            if swarm(i).lfactorial_costs < swarm(i).factorial_costs
                swarm(i).lbvec = swarm(i).uvec;
                swarm(i).lfactorial_costs = swarm(i).factorial_costs;
            end
        end
        
        title(['Generation = ', num2str(generation)]);
        disp(['PSO Generation = ', num2str(generation), ' best factorial costs = ', num2str(bestobj)]);
    end
    if generation >= maxgen
        data_PSO.success = 0;
    end

    data_PSO.wall_clock_time = toc;
    data_PSO.EvBestFitness = EvBestFitness;
    data_PSO.bestInd_data = bestInd_data;
    data_PSO.TotalEvaluations = TotalEvaluations;
    data_PSO.funvalue = funvalue;
    save('data_PSO','data_PSO');
    for i = 1:ntasks
        figure(i) 
        hold all
        plot(EvBestFitness(i,:));
        xlabel('GENERATIONS')
        ylabel(['TASK ', num2str(i), ' OBJECTIVE'])
        legend('PSO')
    end
end