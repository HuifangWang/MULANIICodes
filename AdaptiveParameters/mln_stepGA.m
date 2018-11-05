function [nextPopulation,state] = mln_stepGA(thisScore,thisPopulation,options,state,GenomeLength)
% updated by Huifang Wang, get the next population, 2015/05/28 

%STEPGA Moves the genetic algorithm forward by one generation
%   This function is private to GA.

%   Copyright 2003-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/08/29 08:27:20 $

% how many crossover offspring will there be from each source?
nEliteKids = options.EliteCount;
nXoverKids = round(options.CrossoverFraction * (size(thisPopulation,1) - nEliteKids));
nMutateKids = size(thisPopulation,1) - nEliteKids - nXoverKids;
% how many parents will we need to complete the population?
nParents = 2 * nXoverKids + nMutateKids;

% decide who will contribute to the next generation

% fitness scaling
state.Expectation = mln_fitscalingrank(thisScore,nParents);

% selection. parents are indices into thispopulation
parents = mln_selectionstochunif(state.Expectation,nParents);

% shuffle to prevent locality effects. It is not the responsibility
% if the selection function to return parents in a "good" order so
% we make sure there is a random order here.
parents = parents(randperm(length(parents)));

[~,k] = sort(thisScore,'descend');

% Everyones parents are stored here for genealogy display
state.Selection = [k(1:options.EliteCount),parents];

% here we make all of the members of the next generation
eliteKids  = thisPopulation(k(1:options.EliteCount),:);
xoverKids  = mln_crossoverscattered(parents(1:2*nXoverKids),GenomeLength,thisPopulation);
mutateKids = mln_GA_multation(parents((1 + 2 * nXoverKids):end),options.MutationGens,GenomeLength,thisScore,thisPopulation,[options.PopInitRangeMin;options.PopInitRangeMax]);

% group them into the next generation
nextPopulation = [ eliteKids ; xoverKids ; mutateKids ];


function expectation = mln_fitscalingrank(scores,nParents)
% FITSCALINGRANK Rank based fitness scaling (single objective only).
%   EXPECTATION = FITSCALINGRANK(SCORES,NPARENTS) calculates the
%   EXPECTATION using the SCORES and number of parents NPARENTS.
%   This relationship can be linear or nonlinear.
%
%   Example:
%   Create an options structure using FITSCALINGRANK as 
%   the fitness scaling function
%     options = gaoptimset('FitnessScalingFcn',@fitscalingrank);

%   Copyright 2003-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2010/02/08 22:34:16 $

scores = -scores(:);
[~,i] = sort(scores);

expectation = zeros(size(scores));
expectation(i) = 1 ./ ((1:length(scores))  .^ 0.5);

expectation = nParents * expectation ./ sum(expectation);

function parents = mln_selectionstochunif(expectation,nParents)
%SELECTIONSTOCHUNIF Choose parents using stochastic universal sampling (SUS).
%   PARENTS = SELECTIONSTOCHUNIF(EXPECTATION,NPARENTS,OPTIONS) chooses the 
%   PARENTS using roulette wheel and uniform sampling, based on EXPECTATION 
%   and number of parents NPARENTS. 
%
%   Example:
%   Create an options structure using SELECTIONSTOCHUNIF as the selection
%   function
%     options = gaoptimset('SelectionFcn', @selectionstochunif);

%   Given a roulette wheel with a slot for each expectation whose size is 
%   equal to the expectation. We then step through the wheel in equal size 
%   steps, so as to cover the entire wheel in nParents steps. At each step, 
%   we create a parent from the slot we have landed in. This mechanism is 
%   fast and accurate.

%   Copyright 2003-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/08/29 08:25:36 $

expectation = expectation(:,1);
wheel = cumsum(expectation) / nParents;

parents = zeros(1,nParents);

% we will step through the wheel in even steps.
stepSize = 1/nParents;

% we will start at a random position less that one full step
position = rand * stepSize;

% a speed optimization. Position is monotonically rising.
lowest = 1; 

for i = 1:nParents % for each parent needed,
    for j = lowest:length(wheel) % find the wheel position
        if(position < wheel(j)) % that this step falls in.
            parents(i) = j;
            lowest = j;
            break;
        end
    end
    position = position + stepSize; % take the next step.
end


function xoverKids  = mln_crossoverscattered(parents,GenomeLength,thisPopulation)
%CROSSOVERSCATTERED Position independent crossover function.
%   XOVERKIDS = CROSSOVERSCATTERED(PARENTS,OPTIONS,GENOMELENGTH, ...
%   FITNESSFCN,SCORES,THISPOPULATION) creates the crossover children XOVERKIDS
%   of the given population THISPOPULATION using the available PARENTS.
%   In single or double point crossover, genomes that are near each other tend
%   to survive together, whereas genomes that are far apart tend to be
%   separated. The technique used here eliminates that effect. Each gene has an
%   equal chance of coming from either parent. This is sometimes called uniform
%   or random crossover.
%
%   Example:
%    Create an options structure using CROSSOVERSCATTERED as the crossover
%    function
%     options = gaoptimset('CrossoverFcn' ,@crossoverscattered);

%   Copyright 2003-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/08/29 08:24:20 $


% How many children to produce?
nKids = floor(length(parents)/2);
% Extract information about linear constraints, if any
%linCon = options.LinearConstr;
%constr = ~isequal(linCon.type,'unconstrained');
% Allocate space for the kids
xoverKids = zeros(nKids,GenomeLength);

% To move through the parents twice as fast as thekids are
% being produced, a separate index for the parents is needed
index = 1;
% for each kid...
for i=1:nKids
    % get parents
    r1 = parents(index);
    index = index + 1;
    r2 = parents(index);
    index = index + 1;
    % Randomly select half of the genes from each parent
    % This loop may seem like brute force, but it is twice as fast as the
    % vectorized version, because it does no allocation.
    for j = 1:GenomeLength
        if(rand > 0.5)
            xoverKids(i,j) = thisPopulation(r1,j);
        else
            xoverKids(i,j) = thisPopulation(r2,j);
        end
    end
    % Make sure that offspring are feasible w.r.t. linear constraints
%     if constr
%         feasible  = isTrialFeasible(xoverKids(i,:)',linCon.Aineq,linCon.bineq,linCon.Aeq, ...
%             linCon.beq,linCon.lb,linCon.ub,sqrt(options.TolCon));
%         if ~feasible % Kid is not feasible
%             % Children are arithmetic mean of two parents (feasible w.r.t
%             % linear constraints)
%             alpha = rand;
%             xoverKids(i,:) = alpha*thisPopulation(r1,:) + ...
%                 (1-alpha)*thisPopulation(r2,:);
%         end
%     end
end


function mutateKids = mln_GA_multation(parents, MutationGens,GenomeLength,thisScore,thisPopulation,ranges)

Nmut=length(parents);
mutateKids=thisPopulation(parents,:);

for imut = 1: Nmut
    if thisScore(parents(imut))<0.7
        NMutationGens=floor(GenomeLength/2);
    else
        NMutationGens=MutationGens;
    end
    
    indMutate=randperm(GenomeLength,NMutationGens);
    for igene=1:NMutationGens
        mutateKids (imut,indMutate(igene)) = rand*(abs(diff(ranges(:,indMutate(igene)))))+ranges(1,indMutate(igene));
    end
end
        


   
