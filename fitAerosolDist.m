% Function to fit aerosol distirbution

function [A, mu, sigma] = fitAerosolDist(diameters, densities, varargin)

    p = inputParser;
    addOptional(p,'initPop',[]);
    addOptional(p,'fitType','counts');
    parse(p,varargin{:});

    initPop_d = p.Results.initPop;
    initPop_r = log(exp(initPop_d)/2);
    
    popSize = 50;
    if (size(initPop_r,1) > popSize/2)
        initPop_r = initPop_r(end-popSize/2:end,:);
    end
    
    
    fitType = p.Results.fitType;

    % Using models based around radii fitting
    % (http://eodg.atm.ox.ac.uk/user/grainger/research/aerosols.pdf)
    radii = diameters/2;
    diameters_av = (diameters(1:6)+diameters(2:7))/2; % Mean diameter in each bin
    binSizes = (diameters(2:end) - diameters(1:end-1));
    counts = densities .* repmat(binSizes,size(densities,1),1); % Since we are integrating the PDF can just use counts directly
    vols = 4/3*pi*(diameters_av/2).^3 * (1e-6)^3;


    % set up a function to numerically integrate over the bin size range
    nIntSteps = 4;
    
    if strcmpi(fitType,'counts')
        intFun = @(y) trapz(lognpdf(linspace(y(1),y(2),nIntSteps),y(3),y(4)))*(y(2)-y(1))/(nIntSteps-1); % For particle count
    elseif strcmpi(fitType,'volume')
        intFun = @(y) trapz(1./linspace(y(1),y(2),nIntSteps).*lognpdf(linspace(y(1),y(2),nIntSteps),y(3),y(4)))*(y(2)-y(1))/(nIntSteps-1); % For total volume
    end
    
    prior = @(x) 1; % Uninformative prior
    
    objFun = @(x) -1*(counts(1)*log(intFun([radii(1),radii(2),x(1),x(2)]))...
                    + counts(2)*log(intFun([radii(2),radii(3),x(1),x(2)]))...
                    + counts(3)*log(intFun([radii(3),radii(4),x(1),x(2)]))...
                    + counts(4)*log(intFun([radii(4),radii(5),x(1),x(2)]))...
                    + counts(5)*log(intFun([radii(5),radii(6),x(1),x(2)]))...
                    + counts(6)*log(intFun([radii(6),radii(7),x(1),x(2)]))...
                    + log(prior(x)));

    % First run a genetic search to avoid local minima        
    options_ga = optimoptions('ga', 'TolFun',1e-14, 'InitialPopulationMatrix',initPop_r, 'MaxGenerations',1000, 'PopulationSize', popSize, 'Display', 'off');
    bestVal = ga(objFun,2,[],[],[],[],[],[],[],options_ga);
            
   
    % Refine using a local search
    options_fminunc = optimoptions('fminunc','MaxFunctionEvaluations',1e4, 'StepTolerance',1e-10, 'OptimalityTolerance',1e-12, 'Display', 'off');
    bestVal = fminunc(objFun,bestVal,options_fminunc);
    
   
    % Convert back to diameters
    A = sum(counts);
    mu_r = bestVal(1);
    mu = log(exp(mu_r)*2);
    sigma_r = bestVal(2);
    sigma = log(exp(sigma_r)*2);
    
    normConst = logncdf(diameters(7),mu,sigma) - logncdf(diameters(1),mu,sigma);

    plotFit = true;
    
    plotDiams = exp(linspace(-2,3,100));
    plotDiams_av = 0.5*(plotDiams(1:end-1) + plotDiams(2:end));
    plotBinSizes = (plotDiams(2:end) - plotDiams(1:end-1));
  
    if plotFit
        plot(log(diameters_av),densities/A*normConst);
        hold on;
        %plot(log(diameters_av),(logncdf(diameters(2:7),mu,sigma) - logncdf(diameters(1:6),mu,sigma))./binSizes,'r-');
        plot(log(plotDiams_av),(logncdf(plotDiams(2:end),mu,sigma) - logncdf(plotDiams(1:end-1),mu,sigma))./plotBinSizes,'r-');
        hold off;
        pause(0.01)
    end


end