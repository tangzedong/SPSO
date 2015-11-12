function [data_MFPSO, ppt] = MFPSO(Tasks, option)
    clc    
    tic   
%%%初始化参数
    data_MFPSO.success = 1;
    c1 = option.c1;
    c2 = option.c2;
    c3 = option.c3;
    momentum = option.momentum;
    pop = option.pop;
    maxgen = option.maxgen;
    rmp = option.rmp;
    p_il = option.p_il;
    rho = option.alpha;
    mu = 10;
    impulse = 0;
    rmps = [0.2,0.4,0.6,0.8,1];
    if mod(pop,2) ~= 0  %保证种群数是偶数
        pop = pop + 1;
    end
    ntasks = length(Tasks);  %任务数
    if ntasks <= 1
        error('At least 2 tasks required for MFPSO');
    end
    Dim = zeros(1,ntasks);  %每个任务的维数
    for i = 1:ntasks
        Dim(i) = Tasks(i).dims;
    end
    dim = max(Dim);  %取任务中最大的维数作为个体向量的维数
    local_search_opt = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton');   %个体实施局部优化
    
    fnceval_calls = 0;      %计算目标函数执行次数
    calls_per_individual = zeros(1,pop);
    EvBestFitness = zeros(ntasks,maxgen);  %记录最好适应度
    TotalEvaluations = zeros(1,maxgen);
    bestobj = inf*(ones(1,ntasks));      %记录找到的最优目标函数值

%%%初始化种群
    for i = 1 : pop
        swarm(i) = Particle();         
        swarm(i) = initialize(swarm(i),dim);
        swarm(i).skill_factor = 0;
    end
    parfor i = 1 : pop
        [swarm(i),calls_per_individual(i)] = evaluate(swarm(i),Tasks,p_il,ntasks,local_search_opt);
    end
    
    fnceval_calls = fnceval_calls + sum(calls_per_individual);        %累加每个个体的目标函数计算次数
    TotalEvaluations(1) = fnceval_calls; %TotalEvaluations记录每次迭代完成时目标函数的执行次数
    
%%%初始化种群的factorial_cost, factorial_ranks
    factorial_cost = zeros(1,pop);
    for i = 1:ntasks
        for j = 1:pop
            factorial_cost(j) = swarm(j).factorial_costs(i);
        end
        [xxx,y]=sort(factorial_cost);
        swarm = swarm(y);
        for j=1:pop
            swarm(j).factorial_ranks(i)=j; 
        end
        bestobj(i) = swarm(1).factorial_costs(i);   %bestobj(i)记录目标函数(i)的最好值
        oldbest(i) = bestobj(i);
        EvBestFitness(i,1) = bestobj(i);              %EvBestFitness(i)记录每一代的(i)最优值
        bestInd_data(i) = swarm(1);              %bestInd_data(i)记录最优个体
        bestvec(i, :) = swarm(1).uvec;
%         bestfitness(i) = swarm(1).scalar_fitness;
    end
    for i=1:pop
        swarm(i).lscalar_fitness = 0;
        [minfactorial, minindex]=min(swarm(i).factorial_ranks);
        x=find(swarm(i).factorial_ranks == minfactorial);
        equivalent_skills=length(x);
        if equivalent_skills > 1
            swarm(i).skill_factor = x(1+round((equivalent_skills-1)*rand(1)));
            tmp = swarm(i).factorial_costs(swarm(i).skill_factor);
            swarm(i).factorial_costs(1:ntasks) = inf;
            swarm(i).factorial_costs(swarm(i).skill_factor) = tmp;
        else
            swarm(i).skill_factor = minindex;
            tmp=swarm(i).factorial_costs(swarm(i).skill_factor);
            swarm(i).factorial_costs(1:ntasks) = inf;
            swarm(i).factorial_costs(swarm(i).skill_factor) = tmp;
        end
        [minfactorial,minindex] = min(swarm(i).factorial_ranks);
        swarm(i).skill_factor = minindex(1);
        swarm(i).scalar_fitness = 1 / minfactorial(1);
        swarm(i).ex_fitness = swarm(i).scalar_fitness;
    end

    ntotalcall = 1;
    generation=0;
    eps = 0.1;
    while generation <= maxgen %bestobj(1) > 0.1 %
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
        generation;
%%to do 
%%pso的个体速度和位置的更新，选择最优个体，记录局部最优位置
        oldskill_factor = [swarm(:).skill_factor];
        u = rand(1,dim);
        cf = zeros(1,dim);
        cf(u<=0.5)=(2*u(u<=0.5)).^(1/(mu+1));
        cf(u>0.5)=(2*(1-u(u>0.5))).^(-1/(mu+1));
        for i = 1 : pop
            Vel = swarm(i).velocity;
            Vel = momentum * Vel ...
                + c2 .* rand() .* bsxfun(@minus, bestvec(round(swarm(i).skill_factor), :), swarm(i).uvec) ...
                + c1 .* rand() .* bsxfun(@minus, swarm(i).lbvec, swarm(i).uvec);

            s2 = ceil(pop * rand());
            while i == s2
                s2 = ceil(pop * rand());
            end    

            r1 = rand();
            if ( swarm(i).skill_factor ~= oldskill_factor(s2)) && r1 < rmp%
                if r1<rmp
                    nexchange = nexchange + 1;
                    swarm(i).success_exchange = 1;
                end
                splitpoint = ceil(21*rand(1,3));%1:10:dim;
                splitpoint = splitpoint(randperm(length(splitpoint)));
                perrand = randperm(dim);
                perrand = [];
                for t = 1:length(splitpoint)
                    perrand = [perrand splitpoint(t):splitpoint(t)+9];
                end
                uvec1 = swarm(s2).uvec;
                piovat = ceil(rand(1,3)*dim);
                piovat = sort(piovat);
                tmpvec = {};
                tmpvec{1} = uvec1(1:piovat(1));
                for t = 2:length(piovat)
                    tmpvec{t} = uvec1(piovat(t-1)+1:piovat(t));
                end
                tmpvec{t+1} = uvec1(piovat(t)+1:end);
                per = randperm(length(tmpvec));
                uvec1 = [];
                for t = 1:length(piovat)+1
                    uvec1 = [uvec1 tmpvec{t}];
                end
                Vel = Vel + c3 .* rand() .* (swarm(s2).lbvec - swarm(i).lbvec);%rand(1,Dim_multitask);%
%                 swarm(i).uvec = (1-cf).*swarm(i).uvec + cf.*uvec1;
                if rand() < 0.5
                    swarm(i).skill_factor = oldskill_factor(s2);
                end
            end

            swarm(i).velocity = Vel;
            
            swarm(i).velocity(Vel < -Velmax) = -Velmax;%*(1-abs(normrnd(0, 0.1)));
            swarm(i).velocity(Vel > Velmax) = Velmax;%*(1-abs(normrnd(0, 0.1)));

            PS = swarm(i).uvec;

            PS = rho * PS + Vel;
            swarm(i).uvec = PS + 0.01*normrnd(0, 0.1,1,dim);
            swarm(i).uvec(PS < 0) = abs(normrnd(0, 0.1,1,nnz(PS < 0)));
            swarm(i).uvec(PS > 1) = 1-abs(normrnd(0, 0.1,1,nnz(PS > 1)));
        end
        
        %%评估更新后的个体
        parfor i = 1 : pop            
            [swarm(i),calls_per_individual(i)] = evaluate(swarm(i),Tasks,p_il,ntasks,local_search_opt);    
        end
        fnceval_calls = fnceval_calls + sum(calls_per_individual);
        TotalEvaluations(generation) = fnceval_calls;
        
        for i = 1:pop
            ppt(:, ntotalcall) = evaluate_TEST(swarm(i),Tasks,p_il,ntasks,local_search_opt);
            ntotalcall = ntotalcall + 1;
        end
        
        factorial_cost = zeros(1,pop);
        for i = 1:ntasks
            for j = 1:pop
                factorial_cost(j) = swarm(j).factorial_costs(i);
            end
            [~,y] = sort(factorial_cost);
            swarm = swarm(y);
            externalpop((i - 1) * 10 + 1:(i - 1) * 10 + 10) = swarm(1:10);

            for j = 1:pop
                swarm(j).factorial_ranks(i) = j;
            end
            if swarm(1).factorial_costs(i) <= bestobj(i)
                bestvec(i,:) = swarm(1).uvec;
                bestobj(i) = swarm(1).factorial_costs(i);
                bestInd_data(i) = swarm(1);
                data_MFPSO.bestindx = y(1);
                data_MFPSO.nfc = y(1) + generation * pop;
            end
            EvBestFitness(i,generation) = bestobj(i);            
        end
        
        for i=1:pop
            [minfactorial,minindex] = min(swarm(i).factorial_ranks);
            swarm(i).skill_factor = minindex(1);
            swarm(i).scalar_fitness = 1 / minfactorial(1);
            %记录个体历史最优值
            if swarm(i).scalar_fitness > swarm(i).lscalar_fitness
                swarm(i).lbvec = swarm(i).uvec;
                swarm(i).lscalar_fitness = swarm(i).scalar_fitness;
                swarm(i).lfactorial_ranks = swarm(i).factorial_ranks;
                swarm(i).lfactorial_costs = swarm(i).factorial_costs;
            end
            
            if swarm(i).scalar_fitness < swarm(i).ex_fitness && swarm(i).success_exchange == 1
                swarm(i).success_exchange = 0;
            end
            swarm(i).ex_fitness = swarm(i).scalar_fitness;
        end
        
        nsucexchg = 0;
        for i = 1:length(swarm)
            nsucexchg = nsucexchg + swarm(i).success_exchange;
        end
        oldimpulse = impulse;
        if nexchange ~= 0
            impulse = nsucexchg/nexchange;
        else
            impulse = rmp;
        end
        for i = 1:length(externalpop)
            ppt2(:, i) = evaluate_TEST(externalpop(i),Tasks,p_il,ntasks,local_search_opt);
        end
        ppt2n(1, :) = (ppt2(1, :) - min(ppt2(1,:)))./ (max(ppt2(1, :)) - min(ppt2(1,:)));
        ppt2n(2, :) = (ppt2(2, :) - min(ppt2(2,:)))./ (max(ppt2(2, :)) - min(ppt2(2,:)));
        [idxppt2, Cppt2, ~, Dppt2] = kmeans(ppt2n', 2, 'MaxIter', 10);
        if (nnz(Cppt2(1,:) <= Cppt2(2,:)) == ntasks && ...
                nnz(Cppt2(1, :) < Cppt2(2, :)) > 0) || ...
                (nnz(Cppt2(1,:) >= Cppt2(2,:)) == ntasks && ...
                nnz(Cppt2(1, :) > Cppt2(2, :)) > 0)
            rmp = 1;
        else
            rmp = rmp + sign(impulse-oldimpulse) * normrnd(0, 0.01);
        end
        if nnz(bestobj < oldbestobj) >= 0
            oldbestobj(bestobj < oldbestobj) = bestobj(bestobj < oldbestobj);
        else
            rmp = rmps(ceil(length(rmps)*rand()));
        end
        disp(rmp);
        disp(['impulse = ', num2str(impulse)]);
        
        title(['Generation = ', num2str(generation)]);
        disp(['MFPSO Generation = ', num2str(generation), ' best factorial costs = ', num2str(bestobj)]);
    end
    if generation >= maxgen
        data_MFPSO.success = 0;
    end

    data_MFPSO.wall_clock_time = toc;
    data_MFPSO.EvBestFitness = EvBestFitness;
    data_MFPSO.bestInd_data = bestInd_data;
    data_MFPSO.TotalEvaluations = TotalEvaluations;
    save('data_MFPSO','data_MFPSO');
    for i = 1:ntasks
        figure(i) 
        hold all
        plot(EvBestFitness(i,:));
        xlabel('GENERATIONS')
        ylabel(['TASK ', num2str(i), ' OBJECTIVE'])
        legend('MFPSO')
    end
end