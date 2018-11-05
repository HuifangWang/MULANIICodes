function Population = mln_gacreationuniform(GenomeLength,options)
%GACREATIONUNIFORM Creates the initial population for genetic algorithm.
%   POP = GACREATIONUNIFORM(NVARS,FITNESSFCN,OPTIONS) Creates the
%   initial population that GA will then evolve into a solution.
%
%   Example:
%     options = gaoptimset('CreationFcn',@gacreationuniform);

%   Copyright 2003-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2010/11/01 19:36:46 $
%

if strcmpi(options.PopulationType,'custom')
    error(message('globaloptim:gacreationuniform:unknownPopulationType', options.PopulationType));
end

totalPopulation = sum(options.PopulationSize);
initPopProvided = size(options.InitialPopulation,1);
individualsToCreate = totalPopulation - initPopProvided;

if strcmpi(options.PopulationType,'doubleVector')
    % Initialize Population to be created
    Population = zeros(totalPopulation,GenomeLength);
    % Use initial population provided already
    if initPopProvided > 0
        Population(1:initPopProvided,:) = options.InitialPopulation;
    end
    % Create remaining population
        lowerBound = options.PopInitRangeMin;
        span = options.PopInitRangeMax - lowerBound;
        Population(initPopProvided+1:end,:) = repmat(lowerBound,individualsToCreate,1) + ...
            repmat(span,individualsToCreate,1) .* rand(individualsToCreate,GenomeLength);
  
end

if all(isnan(Population))
    error(message('globaloptim:gacreationuniform:populationIsNaN'));
elseif all(isinf(Population))
    error(message('globaloptim:gacreationuniform:populationIsInf'));
end

