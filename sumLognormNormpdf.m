function [logprob, norm_const] = sumLognormNormpdf(x, mu_ln, sig_ln, mu_n, sig_n, nVals)

    startVal_all = mu_n - 10*sig_n;
    endVal_n = mu_n + 10*sig_n;
    endVal_ln = exp(mu_ln + 4*sig_ln);
    
%     startVal = -10000;
%     endVal_ln = 1e5;


    %nVals = 1e4;

    xVals_ln = linspace(min(startVal_all), max(endVal_ln), nVals);
    dx = xVals_ln(2) - xVals_ln(1);
    xVals_n = startVal_all:dx:endVal_n;


    vals_ln = lognpdf(xVals_ln, mu_ln, sig_ln);
    vals_n = normpdf(xVals_n, mu_n, sig_n);

    vals_sum = conv(vals_ln, vals_n, 'same') * dx;
    %vals_sum = vals_ln;

    norm_const = dx*sum(vals_sum);

    %vals_sum = vals_sum/norm_const;
    
    prob = interp1(xVals_ln, vals_sum, x);
    prob(isnan(prob)) = 1e-50; %deal with values outside range
    prob(prob == 0) = 1e-50; %deal with values outside range
    
    logprob = log(prob);

%     figure;
%     subplot(1,3,1);
%     plot(xVals, vals_ln);
% 
%     subplot(1,3,2);
%     plot(xVals, vals_n);
% 
%     subplot(1,3,3);
%     plot(xVals, vals_sum);

end
