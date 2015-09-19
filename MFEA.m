clear
MaxIter = 100;
nDim = 30;
nPop = 200;
alpha = 0.1;
c2 = 1.2;
c1 = 0.9;
momentum = 0.2;

M=orth(randn(nDim,nDim));
Lb=-0.5*ones(1,nDim);
Ub=0.5*ones(1,nDim);
PS = rand(nPop, nDim); 
Vel = 2 * rand(nPop, nDim) - 1;
PStmp = bsxfun(@plus, bsxfun(@times, PS, (Ub - Lb)), Lb);
fitness = evaluate(PStmp, M);
localBests = PS;
[Neg, besti] = min(fitness);
gobalBfit = Neg(1);
gobalBest = PS(besti(1), :);
localBfit = fitness;

i = 0;
while i < MaxIter
    
    Vel = momentum * Vel + c2 * rand() * bsxfun(@minus, gobalBest, PS) + c1 * rand() * bsxfun(@minus, localBests, PS);
    Vel(Vel <= -1) = -1;
    Vel(Vel >= 1) = 1;
    PS = alpha * PS + (MaxIter - i)/MaxIter * Vel;
    PS(PS < 0) = 0;
    PS(PS > 1) = 1;
    PStmp = bsxfun(@plus, bsxfun(@times, PS, (Ub - Lb)), Lb);
    fitness = evaluate(PStmp, M);
    localBests(fitness < localBfit, :) = PS(fitness < localBfit, :);
    localBfit(fitness < localBfit) = fitness(fitness < localBfit);
    
    [Neg, besti] = min(fitness);
    if gobalBfit > Neg(1)
        gobalBest = PS(besti(1), :);
        gobalBfit = Neg(1);
    end
    [Neg] = min(localBfit);
    fittrace(i + 1) = Neg(1);
    i = i + 1;
end

plot(fittrace, 'r-*');