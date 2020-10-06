% Function to fit aerosol distirbution

function [A, mu, sigma, likelihood] = fitAerosolDist(diameters, densities, varargin)

    p = inputParser;
    addOptional(p,'initPop',[]);
    addOptional(p,'fitType','counts');
    parse(p,varargin{:});

    initPop_d = p.Results.initPop;
    
    if (size(initPop_d,1) > 0)
        initPop_r = [log(exp(initPop_d(:,1))/2),initPop_d(:,2)];
    else
        initPop_r = initPop_d;
    end
    
    popSize = 100;
    if (size(initPop_r,1) > popSize/2)
        initPop_r = initPop_r(end-popSize/2:end,:);
    end
    
    
    fitType = p.Results.fitType;

    logDensity = true;
    
    % Using models based around radii fitting
    % (http://eodg.atm.ox.ac.uk/user/grainger/research/aerosols.pdf)
    log_diameters = log(diameters);
    radii = diameters/2;
    log_radii = log(radii);
    
    diameters_av = (diameters(1:6)+diameters(2:7))/2; % Mean diameter in each bin
    log_diameters_av = (log_diameters(1:6)+log_diameters(2:7))/2;
    bin_sizes = (diameters(2:end) - diameters(1:end-1));
    log_bin_sizes = log(diameters(2:7)) - log(diameters(1:6));
    
    if logDensity
        counts = densities .* repmat(log_bin_sizes,size(densities,1),1);
    else
        counts = densities .* repmat(bin_sizes,size(densities,1),1); % Since we are integrating the PDF can just use counts directly
    end
    
    vols = 4/3*pi*(diameters_av/2).^3 * (1e-6)^3;


    % set up a function to numerically integrate over the bin size range
    nIntSteps = 4;
    
    if strcmpi(fitType,'counts')
        if logDensity
            intFun = @(y) trapz(normpdf(linspace(y(1),y(2),nIntSteps),y(3),y(4)))*(y(2)-y(1))/(nIntSteps-1); % For particle count
        else
            intFun = @(y) trapz(lognpdf(linspace(y(1),y(2),nIntSteps),y(3),y(4)))*(y(2)-y(1))/(nIntSteps-1); % For particle count
        end
    elseif strcmpi(fitType,'volume')
        if logDensity
        else
            intFun = @(y) trapz(1./linspace(y(1),y(2),nIntSteps).*lognpdf(linspace(y(1),y(2),nIntSteps),y(3),y(4)))*(y(2)-y(1))/(nIntSteps-1); % For total volume
        end
    end
    
    mu_UB = log(0.3/2); %D = 0.2, sig = 0.6 is a good fit
    sig_UB = 0.9;
    
    mu_LB = log(0.05/2);
    sig_LB = 0.01;
    
    log_prior = @(x) 0; % Uninformative prior
    %log_prior = @(x) log(normpdf(x(1),log(0.36/2),0.1)); % Prior is on the radius...
    log_prior = @(x) log(unifpdf(x(1),mu_LB, mu_UB)) + log(unifpdf(x(2),sig_LB,sig_UB)); % Prior is on the radius...
    
   
    
    if logDensity
%         objFun = @(x) -1*(counts(1)*log(intFun([log_radii(1),log_radii(2),x(1),x(2)]))...
%                         + counts(2)*log(intFun([log_radii(2),log_radii(3),x(1),x(2)]))...
%                         + counts(3)*log(intFun([log_radii(3),log_radii(4),x(1),x(2)]))...
%                         + counts(4)*log(intFun([log_radii(4),log_radii(5),x(1),x(2)]))...
%                         + counts(5)*log(intFun([log_radii(5),log_radii(6),x(1),x(2)]))...
%                         + counts(6)*log(intFun([log_radii(6),log_radii(7),x(1),x(2)]))...
%                         + log_prior(x));
                    
        objFun = @(x) -1* (normApproxMNpdf(x,log_radii,counts,intFun))...
                        - log_prior(x);
                    
        objFun_dir = @(x) log(0*(densities(1) - x(3)*(normcdf(log_radii(2),x(1),x(2)) - normcdf(log_radii(1),x(1),x(2))))^2 ...
                         +(densities(2) - x(3)*(normcdf(log_radii(3),x(1),x(2)) - normcdf(log_radii(2),x(1),x(2))))^2 ...
                         +(densities(3) - x(3)*(normcdf(log_radii(4),x(1),x(2)) - normcdf(log_radii(3),x(1),x(2))))^2 ...
                         +(densities(4) - x(3)*(normcdf(log_radii(5),x(1),x(2)) - normcdf(log_radii(4),x(1),x(2))))^2 ...
                         +(densities(5) - x(3)*(normcdf(log_radii(6),x(1),x(2)) - normcdf(log_radii(5),x(1),x(2))))^2 ...
                         +(densities(6) - x(3)*(normcdf(log_radii(7),x(1),x(2)) - normcdf(log_radii(6),x(1),x(2))))^2);
                    
        log_likelihood = @(x) -1 * (normApproxMNpdf(x,log_radii,counts,intFun));
    else
        objFun = @(x) -1*(counts(1)*log(intFun([radii(1),radii(2),x(1),x(2)]))...
                        + counts(2)*log(intFun([radii(2),radii(3),x(1),x(2)]))...
                        + counts(3)*log(intFun([radii(3),radii(4),x(1),x(2)]))...
                        + counts(4)*log(intFun([radii(4),radii(5),x(1),x(2)]))...
                        + counts(5)*log(intFun([radii(5),radii(6),x(1),x(2)]))...
                        + counts(6)*log(intFun([radii(6),radii(7),x(1),x(2)]))...
                        + log_prior(x));
    end
    
%     mu_test = linspace(-2.0,3,100);
%     sigma_test = linspace(0.01,2,100);
%     
%     for muIdx = 1:size(mu_test,2)
%         for sigmaIdx = 1:size(sigma_test,2)
%             test(muIdx,sigmaIdx) = objFun([mu_test(muIdx),sigma_test(sigmaIdx)]);
%         end
%     end
%     imagesc(mu_test,sigma_test,-1*log(test),[-16,-14]);
%     xlabel('mu');
%     ylabel('sigma');
    
    
    w = warning;
    warning('off','all');
    % First run a genetic search to avoid local minima        
    options_ga = optimoptions('ga', 'TolFun',1e-14, 'InitialPopulationMatrix',initPop_r, 'MaxGenerations',1000, 'PopulationSize', popSize, 'Display', 'off');
    bestVal = ga(objFun,2,[],[],[],[],[mu_LB,sig_LB],[mu_UB,sig_UB],[],options_ga);
    %bestVal = ga(objFun_dir,3,[],[],[],[],[],[],[],options_ga);
   
    % Refine using a local search
    options_fminunc = optimoptions('fminunc','MaxFunctionEvaluations',1e4, 'StepTolerance',1e-10, 'OptimalityTolerance',1e-12, 'Display', 'off');
    bestVal = fminunc(objFun,bestVal,options_fminunc);
    %bestVal = fminunc(objFun_dir,bestVal,options_fminunc);
    warning(w);
    
    likelihood = -1*log_likelihood(bestVal);
   
    % Convert back to diameters
    A = sum(counts);
    mu_r = bestVal(1);
    mu = log(exp(mu_r)*2);
    sigma_r = bestVal(2);
    sigma = sigma_r; % In log space doesn't change as distribution is just shifted
    
    if logDensity
        normConst = normcdf(log_diameters(7),mu,sigma) - normcdf(log_diameters(1),mu,sigma);
    else
        normConst = logncdf(diameters(7),mu,sigma) - logncdf(diameters(1),mu,sigma);
    end

    plotFit = true;
    
    plotDiams = exp(linspace(-0.5,3,300));
    log_plotDiams = log(plotDiams);
    plotDiams_av = 0.5*(plotDiams(1:end-1) + plotDiams(2:end));
    log_plotDiams_av = (log_plotDiams(1:end-1)+log_plotDiams(2:end))/2;
    plotBinSizes = (plotDiams(2:end) - plotDiams(1:end-1));
    log_plotBinSizes = log(plotDiams(2:end)) - log(plotDiams(1:end-1));
  
    if plotFit
        if logDensity
            plot(log_diameters_av,densities/A*normConst);
        else
            plot(log(diameters_av),densities/A*normConst);
        end
        hold on;
        %plot(log(diameters_av),(logncdf(diameters(2:7),mu,sigma) - logncdf(diameters(1:6),mu,sigma))./binSizes,'r-');
        if logDensity
            plot(log_plotDiams_av,(normcdf(log_plotDiams(2:end),mu,sigma) - normcdf(log_plotDiams(1:end-1),mu,sigma))./log_plotBinSizes,'r-');
        else
            plot(log(plotDiams_av),(logncdf(plotDiams(2:end),mu,sigma) - logncdf(plotDiams(1:end-1),mu,sigma))./plotBinSizes,'r-');
        end
        
        %err = (densities/A*normConst) - ((logncdf(diameters(2:end),mu,sigma) - logncdf(diameters(1:end-1),mu,sigma))./bin_sizes);
        %err = abs(err)./(densities/A*normConst);
        %plot(log(diameters_av),err,'k:');
        
        disp(['Log likelihood = ', num2str(log_likelihood(bestVal))]);

        hold off;
        pause(0.01)
    end


end

function prob = normApproxMNpdf(x,log_radii, counts, intFun)

    % https://stats.stackexchange.com/questions/34547/what-is-the-normal-approximation-of-the-multinomial-distribution?noredirect=1&lq=1
    % http://www.stat.umn.edu/geyer/5102/notes/brand.pdf
    probs = [intFun([log_radii(1),log_radii(2),x(1),x(2)]),...
             intFun([log_radii(2),log_radii(3),x(1),x(2)]),...
             intFun([log_radii(3),log_radii(4),x(1),x(2)]),...
             intFun([log_radii(4),log_radii(5),x(1),x(2)]),...
             intFun([log_radii(5),log_radii(6),x(1),x(2)]),...
             intFun([log_radii(6),log_radii(7),x(1),x(2)])];
    
    % Reduce dims.     
%     probs = [probs,1-sum(probs)];
%     
%     P = diag(probs);
%     M = P - probs'*probs;
%     
%     %counts = log(counts);
%     n = sum(counts);
%     
%     mu = n*probs;
%     cov = n*M;
%     
%     cov = cov(1:end-1,1:end-1);
%     mu = mu(1:end-1);
%     
%     prob = -1*0.5*(counts - mu)*inv(cov)*(counts-mu)';% - 0.5*log(abs(det(cov))*(2*pi)^size(counts,2));

    % Ues only some dims
    probs = probs/sum(probs);
    
    P = diag(probs);
    M = P - probs'*probs;
    
    n = sum(counts);
    
    mu = n*probs;
    cov = n*M;
    
    %cov = cov(1:end-1,1:end-1);
    %mu = mu(1:end-1);
    %counts = counts(1:end-1);
    
    prob = -1*0.5*(counts - mu)*inv(cov)*(counts-mu)';% - 0.5*log(abs(det(cov))*(2*pi)^size(counts,2));
   
end