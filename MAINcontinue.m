% This MATLAB R2014b code is for EVOLUTIONARY MULTITASKING across minimization problems. 
% For maximization problems, multiply objective function by -1.

% For suggestions please contact: Abhishek Gupta (Email: abhishekg@ntu.edu.sg or
% agup839@aucklanduni.ac.nz or abhi.nitr2010@gmail.com)

clear
n = 30;
M=orth(randn(n,n));
Load

% LoadTest;
% i = 1;
% %close all
% %% Example 1 - (40-D Rastrigin, 30-D Ackley)
% % % Rastrigin function definition
% % n=3;
% M=orth(randn(n,n));
% % Tasks(1) = testfun(3);
% Tasks(1).dims=n;
% 
% Tasks(1).fncLS=@(fnc, x0, option)fminunc(fnc, x0, option);
% Tasks(1).fnc=@(x)Ackley(x,M);
% Tasks(1).Lb=-50*ones(1,n);
% Tasks(1).Ub=50*ones(1,n);
% 
% % Ackley function definition
% n=20;
% M=eye(n,n);
% Tasks(2).dims=n;
% Tasks(2).fncLS=@(fnc, x0, option)fminunc(fnc, x0, option);
% Tasks(2).fnc=@(x)Rastrigin(x,M);
% Tasks(2).Lb=-50*ones(1,n);
% Tasks(2).Ub=50*ones(1,n);

% Example 2 - (50-D Sphere, 30-D Weierstrass)
% Sphere function definition
% n=30;
% Tasks(1).dims=n;
% M=orth(randn(n,n));
% Tasks(1).fnc=@(x)Sphere(x,M);
% Tasks(1).Lb=-50*ones(1,n);
% Tasks(1).Ub=50*ones(1,n);
% % Rastrigin function definition
% n=30;
% Tasks(2).dims=n;
% M=orth(randn(n,n));
% Tasks(2).fnc=@(x)Ackley(x,M);
% Tasks(2).Lb=-50*ones(1,n);
% Tasks(2).Ub=50*ones(1,n);

%% Example 3 - (40-D Rastrigin, 50-D Ackley, 20-D Sphere)
% % Rastrigin function definition
% n=40;
% Tasks(1).dims=n;
% M=orth(randn(n,n));
% Tasks(1).fnc=@(x)Rastrigin(x,M);
% Tasks(1).Lb=-5*ones(1,n);
% Tasks(1).Ub=5*ones(1,n);
% % Ackley function definition
% n=50;
% Tasks(2).dims=n;
% M=orth(randn(n,n));
% Tasks(2).fnc=@(x)Ackley(x,M);
% Tasks(2).Lb=-32*ones(1,n);
% Tasks(2).Ub=32*ones(1,n);
% % Sphere function definition
% n=20;
% Tasks(3).dims=n;
% M=eye(n,n);
% Tasks(3).fnc=@(x)Sphere(x,M);
% Tasks(3).Lb=-100*ones(1,n);
% Tasks(3).Ub=100*ones(1,n);

%% Example 5 - Mirror Functions - (30-D Rastrigin, 30-D Rastrigin)
% % Rastrigin function definition
% n=30;
% Tasks(1).dims=n;
% M=orth(randn(n,n));
% Tasks(1).fnc=@(x)Rastrigin(x,M);
% Tasks(1).Lb=-5*ones(1,n);
% Tasks(1).Ub=5*ones(1,n);
% Tasks(2) = Tasks(1);

%% Example 6 - Mirror Functions - (30-D Ackley, 30-D Ackley)
% % Ackley function definition
% n=30;
% Tasks(1).dims=n;
% M=orth(randn(n,n));
% Tasks(1).fnc=@(x)Ackley(x,M);
% Tasks(1).Lb=-32*ones(1,n);
% Tasks(1).Ub=32*ones(1,n);
% Tasks(2)=Tasks(1);

%% Calling the solvers
% For large population sizes, consider using the Parallel Computing Toolbox
% of MATLAB.
% Else, program can be slow.
option.pop=100; % population size
option.maxgen=200; % generation count
option.p_il = 1; % probability of individual learning (BFGA quasi-Newton Algorithm) --> Indiviudal Learning is an IMPORTANT component of the MFEA.
option.rmp = 0.2; % random mating probability
option.c1 = 0.8;
option.c2 = 1.2;
option.c3 = 0.8;
% option.c1 = 2.3;
% option.c2 = 1.3;
% option.c3 = 1.8;
option.alpha = 1;
option.momentum = 0.5;
% for t = 1:length(testfun)
%     Tasks(i).dims=n;
%     Tasks(i) = testfun(i);
%     
%     for i = 1:2
%         [data_MFPSO]=MFPSO(Tasks,option);
%         nfc(i,t) = data_MFPSO.nfc;
%         success(i,t) = data_MFPSO.success;
%     end
% end
% save('experience1mfo.mat','nfc','success');
% [data_PSO]=PSO(Tasks(1),option);
for i = 1:5
    sumnfc = 0;
    sumsuc = 0;
    clear resultdata
    Tasks(1) = testfun(i);
    for t = 1:5
        [data_PSO]=PSO(Tasks(1),option);
        data_PSO.name1 = Tasks(1).funname;
        resultdata(t) = data_PSO;
    end
    result{i} = resultdata;
end
save('data.mat', 'result');
pop = option.pop; % population size
gen = option.maxgen; % generation count
p_il = option.p_il; % probability of individual learning (BFGA quasi-Newton Algorithm) --> Indiviudal Learning is an IMPORTANT component of the MFEA.
selection_pressure = 'roulette wheel'; % choose either 'elitist' or 'roulette wheel'
rmp = option.rmp; % random mating probability
% "task_for_comparison_with_SOO" compares performance of corresponding task in MFO with SOO.
% For Instance, In EXAMPLE 1 ...
% "task_for_comparison_with_SOO" = 1 --> compares 40-D Rastrin in MFO with 40-D
% Rastrigin in SOO.
% "task_for_comparison_with_SOO" = 2 --> compares 30D Ackley in MFO with
% 30D Ackley in SOO.
% task_for_comparison_with_SOO = 2;
% data_SOO=SOO(Tasks(task_for_comparison_with_SOO),task_for_comparison_with_SOO,pop,gen,selection_pressure,p_il);     
