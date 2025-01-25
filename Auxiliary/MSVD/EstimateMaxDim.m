function lMaxDim = EstimateMaxDim( lDeltas, lLargeScales,lMeanS,lMeanGapS )

epsnoise = 1e-3;

lLargeGapIdxs   = almostmax( lMeanGapS(lLargeScales,:), 0.5 );
lLargeGapIdxs   = unique( lLargeGapIdxs(:,1:min([size(lLargeGapIdxs,2),3])) );
lLargeGapSorted = sort(lLargeGapIdxs,'descend');
lRoughYScale    = max(max(lMeanS));
lRoughXScale    = max(lDeltas);
lMaxDim         = 0;

% Look for a linearly increasing singular value with a large gap, across large scales
for k = 1:length(lLargeGapSorted),
    lLinFit = polyfit(lDeltas(lLargeScales),lMeanS(lLargeScales,lLargeGapSorted(k)),1);
    if lLinFit(1)>epsnoise*lRoughYScale/lRoughXScale,
        lMaxDim = lLargeGapSorted(k);
        break;
    end;
end;

lMaxDim = max([max([max(max(lLargeGapSorted)),min([lMaxDim,size(lMeanS,2)-1])]),1]);

return;