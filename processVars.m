function [varsProcessed, varsProcessedCats] = processVars(var_raw, targetCatNames, varargin)

    nCats = size(targetCatNames,2);
    
    if nargin == nCats + 2
        var_temp = zeros(size(var_raw));
        for catIdx = 1:nCats
            currentList = varargin{catIdx};
            
            tempVals = false(size(var_raw));
            for stringIdx = 1:size(currentList,2)
                currentString = currentList{stringIdx};
                
                if strcmpi(currentString, '')
                    tempVals = tempVals | cellfun(@isempty,var_raw);
                    emptyCat = catIdx;
                else
                    temp_raw = regexpi(var_raw, currentString);
                    tempVals = tempVals | ~cellfun(@isempty,temp_raw);
                end
            end
            
            var_temp(tempVals) = catIdx;
        end

        if any(var_temp == 0)
            warning('check categorisations!');
        end
        var_temp(var_temp == 0) = emptyCat;
        varsProcessed = categorical(var_temp,1:nCats,targetCatNames);
        varsProcessedCats = categories(varsProcessed);
    else
        error('Not enough category labels!');
    end
    
end