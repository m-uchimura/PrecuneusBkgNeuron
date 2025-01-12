function [pvals,clustIdx]=calcP_clu(ZI,I_SHUFF,nIter)
if 1%doCluster
    thresh = norminv(.95);
    h0 = ZI>=thresh;
    [clustIdx, ~, clustMass] = getClust(h0,ZI);
    nullMass = nan(nIter,1);
    for iIter = 1:nIter
        tmpPhaseMI = I_SHUFF(iIter,:);
        tmpNull = I_SHUFF;
        tmpNull(iIter,:) = [];
        tmpNull = [ZI'; tmpNull]; %#ok<AGROW>
        tmpPhaseMI = ( tmpPhaseMI-mean(tmpNull,1) ) ./ std(tmpNull,[],1);
       
        tmpH0 = abs(tmpPhaseMI)>=thresh;
        [~, ~, tmpMass] = getClust(tmpH0,tmpPhaseMI);
        if isempty(tmpMass)
            nullMass(iIter) = 0;
        else
            nullMass(iIter) = max(abs(tmpMass));
        end
    end
    
    pvals = sum(nullMass>=abs(clustMass)')./nIter;
else pvals=NaN; clustIdx=NaN;
end