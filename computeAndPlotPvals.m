function [pMean, pSig] = computeAndPlotPvals(d1, d2, noiseMean, noiseStd, cat1Idx, cat2Idx, varargin)

    p = inputParser;
    addOptional(p,'pMeanIn',[]);
    addOptional(p,'pSigIn',[]);
    parse(p,varargin{:});

    pMeanIn = p.Results.pMeanIn;
    pSigIn = p.Results.pSigIn;
    
    if isempty(pMeanIn) || isempty(pSigIn)

         reject1 = d1 < -3*noiseStd;
         reject2 = d2 < -3*noiseStd;

         d1 = d1(~reject1);
         d2 = d2(~reject2);

         d1_trunc = d1(d1>0);
         d2_trunc = d2(d2>0);

         if ~isempty(d1_trunc)
            paramEst1 = lognfit(d1_trunc);
         else
            paramEst1 = [0,1];
         end

         if ~isempty(d2_trunc)
            paramEst2 = lognfit(d2_trunc);
         else
            paramEst2 = [0,1];
         end


         objFun = @(x) -1*sum(sumLognormNormpdf(d1, x(1),x(2),0,noiseStd,1e5));
         x0 = paramEst1;
         LB = [0,0];
         UB = [2*paramEst1(1), 2*paramEst1(2)];
         options_fmincon = optimoptions('fmincon','MaxFunctionEvaluations',1e4, 'StepTolerance',1e-12, 'OptimalityTolerance',1e-12, 'Display', 'off');
         bestVal = fmincon(objFun,x0,[],[],[],[], LB , UB,[],options_fmincon);
         paramEst1 = bestVal;

         objFun = @(x) -1*sum(sumLognormNormpdf(d2, x(1),x(2),0,noiseStd,1e5));
         x0 = paramEst2;
         LB = [0,0];
         UB = [2*paramEst2(1), 2*paramEst2(2)];
         options_fmincon = optimoptions('fmincon','MaxFunctionEvaluations',1e4, 'StepTolerance',1e-12, 'OptimalityTolerance',1e-12, 'Display', 'off');
         bestVal = fmincon(objFun,x0,[],[],[],[], LB , UB,[],options_fmincon);
         paramEst2 = bestVal;

         [pMean, pSig] = computeSignificance(d1, d2, noiseMean, noiseStd, paramEst1(1), paramEst1(2), paramEst2(1), paramEst2(2));
    else
        pMean = pMeanIn;
        pSig = pSigIn;
    end

    pMean = min([pMean, 1-pMean]);
    pSig = min([pSig, 1-pSig]);
    
    if (pMean <= 0.05) || (pSig <= 0.05)
        tempYlim = ylim;
        dY = 0.08*tempYlim(2);
        yVal = max([d1(:); d2(:)]) + dY*(cat2Idx - cat1Idx);
        xVal = [cat1Idx, cat2Idx];
        ctr = mean(xVal);

        hold on;
        plot(xVal, [yVal, yVal], '-b', 'LineWidth', 2);
        oldYlim = ylim;
        newYMax = max([yVal*1.1, oldYlim(2)]);
        ylim([oldYlim(1), newYMax]);
        text(ctr, yVal+dY, ['p_{\mu}=',num2str(pMean), ', p_{\sigma}', num2str(pSig)], 'HorizontalAlignment', 'Center', 'Color', 'blue');
        hold off;
    end

         
end