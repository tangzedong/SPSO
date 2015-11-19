function [data_PSO, ppt] = PSO(Tasks, option)
    clc    
    tic   
%%%��ʼ������
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
    if mod(pop,2) ~= 0  %��֤��Ⱥ����ż��
        pop = pop + 1;
    end
    ntasks = 1;
    dim = Tasks(1).dims;  %ȡ����������ά����Ϊ����������ά��
    local_search_opt = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton');   %����ʵʩ�ֲ��Ż�
    
    fnceval_calls = 0;      %����Ŀ�꺯��ִ�д���
    calls_per_individual = zeros(1,pop);
    EvBestFitness = zeros(ntasks,maxgen);  %��¼�����Ӧ��
    TotalEvaluations = zeros(1,maxgen);
    bestobj = inf;      %��¼�ҵ�������Ŀ�꺯��ֵ
    funvalue = -1*ones(ntasks,pop);
%%%��ʼ����Ⱥ
    for i = 1 : pop
        swarm(i) = Particle();         
        swarm(i) = initialize(swarm(i),dim);
    end
    for i = 1 : pop
        [swarm(i),calls_per_individual(i)] = evaluate_SOO(swarm(i),Tasks,p_il,local_search_opt);
        funvalue(ntotalcall) = swarm(i).factorial_costs;
        ntotalcall = ntotalcall + 1;
    end
    
    fnceval_calls = fnceval_calls + sum(calls_per_individual);        %�ۼ�ÿ�������Ŀ�꺯���������
    TotalEvaluations(1) = fnceval_calls; %TotalEvaluations��¼ÿ�ε������ʱĿ�꺯����ִ�д���
    
%%%��ʼ����Ⱥ��factorial_cost, factorial_ranks
    factorial_cost = zeros(1,pop);
    for j = 1:pop
        factorial_cost(j) = swarm(j).factorial_costs;
    end
    [xxx,y]=sort(factorial_cost);
    swarm = swarm(y);
    for j=1:pop
        swarm(j).factorial_ranks(i)=j;
    end
    bestobj = swarm(1).factorial_costs;   %bestobj(i)��¼Ŀ�꺯��(i)�����ֵ
    oldbestobj = bestobj;
    EvBestFitness = bestobj;              %EvBestFitness(i)��¼ÿһ����(i)����ֵ
    bestInd_data = swarm(1);              %bestInd_data(i)��¼���Ÿ���
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
%%pso�ĸ����ٶȺ�λ�õĸ��£�ѡ�����Ÿ��壬��¼�ֲ�����λ��

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
        
        %%�������º�ĸ���
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
            %��¼������ʷ����ֵ
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