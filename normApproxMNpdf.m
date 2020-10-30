
function prob = normApproxMNpdf(x,log_radii, counts, intFun, varargin)

    p = inputParser;
    addOptional(p,'tubeCorrection',[]);
    
    parse(p,varargin{:});

    tubeCorrection = p.Results.tubeCorrection;
    if ~isempty(tubeCorrection)
        hasTubeCorrection = true;
    else
        hasTubeCorrection = false;
    end

    % https://stats.stackexchange.com/questions/34547/what-is-the-normal-approximation-of-the-multinomial-distribution?noredirect=1&lq=1
    % https://stats.stackexchange.com/questions/2397/asymptotic-distribution-of-multinomial
    % http://www.stat.umn.edu/geyer/5102/notes/brand.pdf
    probs = [];
    for k=1:size(log_radii,2)-1
        if hasTubeCorrection
            probs = [probs, tubeCorrection(k)*intFun([log_radii(k),log_radii(k+1),x])];
        else
            probs = [probs, intFun([log_radii(k),log_radii(k+1),x])];
        end
    end
    
%     probs = [intFun([log_radii(1),log_radii(2),x]),...
%              intFun([log_radii(2),log_radii(3),x]),...
%              intFun([log_radii(3),log_radii(4),x]),...
%              intFun([log_radii(4),log_radii(5),x]),...
%              intFun([log_radii(5),log_radii(6),x]),...
%              intFun([log_radii(6),log_radii(7),x])];
    
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
    
    if (any(isnan(probs)))
        prob = -1e30;
        return;
    end
    
    P = diag(probs);
    M = P - probs'*probs;
    
    n = sum(counts);
    
    mu = n*probs;
    cov = n*M;
    
    %cov = cov(1:end-1,1:end-1);
    %mu = mu(1:end-1);
    %counts = counts(1:end-1);
    
    prob_test = -1*0.5*(counts - mu)*inv(cov)*(counts-mu)';
    %normAttempt = 0.5*log(abs(det(cov))*(2*pi)^size(counts,2));
    if (isnan(prob_test))
        prob = -1e40;
        return;
    end
    
    
    % Test
    prob = -1*0.5*(counts - mu)*pinv(cov)*(counts-mu)';
    e = eig(cov);
    e(e < 1e-9 * max(e)) = 1;
    det_est = prod(e);
    normAttempt = 0.5*log(abs(det_est)*(2*pi)^size(counts,2));
    

    prob = prob - normAttempt;
    
 
   
end