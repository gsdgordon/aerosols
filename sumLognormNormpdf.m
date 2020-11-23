function [logprob, norm_const] = sumLognormNormpdf(x, mu_ln, sig_ln, mu_n, sig_n, nVals, maxPossibleVal)

    startVal_all = mu_n - 10*sig_n;
    endVal_n = mu_n + 10*sig_n;
    endVal_ln = exp(mu_ln + 4*sig_ln);
    %endVal_ln = min([endVal_ln, maxPossibleVal]);
    
%     startVal = -10000;
%     endVal_ln = 1e5;

    xVals_ln = linspace(startVal_all, endVal_ln, nVals);
    dx = xVals_ln(2) - xVals_ln(1);
    xVals_n = startVal_all:dx:endVal_n;
    
    vals_ln = lognpdf(xVals_ln, mu_ln, sig_ln);
    vals_n = normpdf(xVals_n, mu_n, sig_n);

% 
%     % Quick compute of vals

%     
%     [xConv, xIn] = meshgrid(xVals_ln, x);
%     temp1 = xIn > xConv;
%     temp2 = xIn < xConv;
%     temp3 = (circshift(temp1,[0,1]) & temp2) | (temp1 & circshift(temp2,[0,-1]));
%     
%     nConv1 = fliplr(vals_n);
%     nConv2 = padarray(vals_ln,[0,size(vals_n,2)-1], 'both');
%     
%     estVal = zeros(size(x,1),1);
%     for k=1:size(x,1)
%           
%         idx1 = find(temp3(k,:),1); 
%         sectionStart = idx1+size(vals_n,2);
%         sectionEnd = sectionStart+size(nConv1,2) - 1;
%         
%         val1 = nConv1 * nConv2(sectionStart:sectionEnd)';
%         val2 = nConv1 * nConv2(sectionStart+1:sectionEnd+1)';
%         
%         estVal(k,1) = interp1([xVals_ln(idx1), xVals_ln(idx1+1)], [val1, val2], x(k))*dx;
%     end


%     
%     %% Estimate norm. constant
%     norm_start = exp(mu_ln - 3*sig_ln) - mu_n - 3*sig_n; %FIX what if mu_n is not equal to zero?
%     if mu_n ~= 0
%         error('Code not valid if mu_n is not equal to zero');
%     end
%     norm_end = exp(mu_ln + 3*sig_ln) + mu_n + 3*sig_n;
%     x_norm = linspace(norm_start, norm_end,200)';
%     
%     [xConv, xIn] = meshgrid(xVals_ln, x_norm);
%     temp1 = xIn > xConv;
%     temp2 = xIn < xConv;
%     temp3 = (circshift(temp1,[0,1]) & temp2) | (temp1 & circshift(temp2,[0,-1]));
%     
%     estVal_norm = zeros(size(x_norm,1),1);
%     for k=1:size(x_norm,1)
%           
%         idx1 = find(temp3(k,:),1); 
%         sectionStart = idx1+size(vals_n,2);
%         sectionEnd = sectionStart+size(nConv1,2) - 1;
%         
%         val1 = nConv1 * nConv2(sectionStart:sectionEnd)';
%         val2 = nConv1 * nConv2(sectionStart+1:sectionEnd+1)';
%         
%         estVal_norm(k,1) = interp1([xVals_ln(idx1), xVals_ln(idx1+1)], [val1, val2], x_norm(k))*dx;
%     end
    
        
        
    vals_sum = conv(vals_ln, vals_n, 'same') * dx;
    %vals_sum = vals_ln;

    norm_const = dx*sum(vals_sum);
    %norm_const = sum(estVal_norm)*dx;

    %vals_sum = vals_sum/norm_const;
    
    prob = interp1(xVals_ln, vals_sum, x);
    %prob = estVal;
    prob(isnan(prob)) = 1e-50; %deal with values outside range
    prob(prob == 0) = 1e-50; %deal with values outside range
    
    logprob = log(prob);
    
    if (any(isinf(logprob)))
        a = 1;
    end

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
