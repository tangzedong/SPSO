function data_MFPSO = MFPSO(Tasks,pop,maxgen,selection_process,rmp,p_il)
    clc    
    tic   
%%%��ʼ������
    if mod(pop,2) ~= 0  %��֤��Ⱥ����ż��
        pop = pop + 1;
    end
    ntasks=length(Tasks);  %������
    if ntasks <= 1
        error('At least 2 tasks required for MFPSO');
    end
    Dim=zeros(1,ntasks);  %ÿ�������ά��
    for i=1:ntasks
        Dim(i)=Tasks(i).dims;
    end
    Dim_multitask=max(Dim);  %ȡ����������ά����Ϊ����������ά��
    local_search_opt = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton');   %����ʵʩ�ֲ��Ż�
    
    fnceval_calls = 0;      %����Ŀ�꺯��ִ�д���
    calls_per_individual=zeros(1,pop);
    EvBestFitness = zeros(ntasks,maxgen);  %��¼�����Ӧ��
    TotalEvaluations=zeros(1,maxgen);
    bestobj=inf*(ones(1,ntasks));      %��¼�ҵ�������Ŀ�꺯��ֵ

%%%��ʼ����Ⱥ
    for i = 1 : pop
        swarm(i) = Particle();         
        swarm(i) = initialize(swarm(i),Dim_multitask);
        swarm(i).skill_factor=0;
    end
    parfor i = 1 : pop
        [swarm(i),calls_per_individual(i)] = evaluate(swarm(i),Tasks,p_il,ntasks,local_search_opt);
    end
    
    fnceval_calls=fnceval_calls + sum(calls_per_individual);        %�ۼ�ÿ�������Ŀ�꺯���������
    TotalEvaluations(1)=fnceval_calls; %TotalEvaluations��¼ÿ�ε������ʱĿ�꺯����ִ�д���
    
%%%��ʼ����Ⱥ��factorial_cost, factorial_ranks
    factorial_cost=zeros(1,pop);
    for i = 1:ntasks
        for j = 1:pop
            factorial_cost(j)=swarm(j).factorial_costs(i);
        end
        [xxx,y]=sort(factorial_cost);
        swarm=swarm(y);
        for j=1:pop
            swarm(j).factorial_ranks(i)=j; 
        end
        bestobj(i)=swarm(1).factorial_costs(i);   %bestobj(i)��¼Ŀ�꺯��(i)�����ֵ
        EvBestFitness(i,1)=bestobj(i);              %EvBestFitness(i)��¼ÿһ����(i)����ֵ
        bestInd_data(i)=swarm(1);              %bestInd_data(i)��¼���Ÿ���
    end
    for i=1:pop
        [minfactorial,minindex]=min(swarm(i).factorial_ranks);
        x=find(swarm(i).factorial_ranks == minfactorial);
        equivalent_skills=length(x);
        if equivalent_skills>1
            swarm(i).skill_factor=x(1+round((equivalent_skills-1)*rand(1)));
            tmp=swarm(i).factorial_costs(swarm(i).skill_factor);
            swarm(i).factorial_costs(1:ntasks)=inf;
            swarm(i).factorial_costs(swarm(i).skill_factor)=tmp;
        else
            swarm(i).skill_factor=minindex;
            tmp=swarm(i).factorial_costs(swarm(i).skill_factor);
            swarm(i).factorial_costs(1:ntasks)=inf;
            swarm(i).factorial_costs(swarm(i).skill_factor)=tmp;
        end
    end
        
    c1 = 0.2; %learning rate 1 learn to historial best particle
    c2 = 1.2; %learning rate 2 learn to global best particle
    generation=1;
    while generation <= maxgen 
        generation = generation + 1;
%%to do 
%%����pso�ĸ����ٶȺ�λ�õĸ��£�ѡ�����Ÿ��壬��¼�ֲ�����λ��
        for i = 1 : pop
            Vel = swarm(i).velocity;
            Vel = momentum * Vel ...
                + c2 * rand() * bsxfun(@minus, gobalBest, swarm(i).uvec) ...
                + c1 * rand() * bsxfun(@minus, localBests, swarm(i).uvec);
            swarm(i).velocity = Vel;
            swarm(i).velocity(Vel <= -1) = -1;
            swarm(i).velocity(Vel >= 1) = 1;
            
            PS = swarm(i).uvec;
            PS = alpha * PS + (MaxIter - i)/MaxIter * Vel;
            swarm(i).uvec = PS;
            swarm(i).uvec(PS < 0) = 0;
            swarm(i).uvec(PS > 1) = 1;
        end
        
        %%���������ɵĸ���
        parfor i = 1 : pop            
            [swarm(i),calls_per_individual(i)] = evaluate(swarm(i),Tasks,p_il,ntasks,localsearchopt);    
        end             
        fnceval_calls=fnceval_calls + sum(calls_per_individual);
        TotalEvaluations(generation)=fnceval_calls;
        
        factorial_cost=zeros(1,pop);
        for i = 1:ntasks
            for j = 1:pop
                factorial_cost(j)=swarm(j).factorial_costs(i);
            end
            [~,y]=sort(factorial_cost);
            swarm=swarm(y);
            for j=1:2*pop
                swarm(j).factorial_ranks(i)=j;
            end
            if swarm(1).factorial_costs(i)<=bestobj(i)
                bestobj(i)=swarm(1).factorial_costs(i);
                bestInd_data(i)=swarm(1);
            end
            EvBestFitness(i,generation)=bestobj(i);            
        end
        for i=1:pop
            [minfactorial,minindex]=min(swarm(i).factorial_ranks);
            swarm(i).skill_factor=minindex;
            swarm(i).scalar_fitness=1/minfactorial;
            if swarm(i).lskill_fitness < swarm(i).skill_fitness
                swarm(i).lbvec = swarm(i).uvec;
                swarm(i).lfactorial_rank = swarm(i).factorial_rank;
                swarm(i).lfactorial_cost = swarm(i).factorial_cost;
            end
            if bestfitness < swarm(i).skill_fitness
                bestvec = swarm(i).uvec;
                bestfitness = swarm(i).skill_fitness;
            end
        end   
        
        disp(['MFPSO Generation = ', num2str(generation), ' best factorial costs = ', num2str(bestobj)]);
    end
    data_MFPSO.wall_clock_time=toc;
    data_MFPSO.EvBestFitness=EvBestFitness;
    data_MFPSO.bestInd_data=bestInd_data;
    data_MFPSO.TotalEvaluations=TotalEvaluations;
    save('data_MFPSO','data_MFPSO');
    for i=1:ntasks
        figure(i)
        hold on
        plot(EvBestFitness(i,:));
        xlabel('GENERATIONS')
        ylabel(['TASK ', num2str(i), ' OBJECTIVE'])
        legend('MFPSO')
    end
end