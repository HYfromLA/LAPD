function almostmaxidxs = almostmax( v, thres )

lSortedGaps = zeros(size(v));

% Sort the gaps at each scale
for k = 1:size(v,1),
    [lSortedValues,lSortedGaps(k,:)] = sort(v(k,:),'descend');
    lSmallIdxs = find(lSortedValues<=thres*lSortedValues(1));
    lSortedGaps(k,lSmallIdxs) = lSortedGaps(k,1);
end;

almostmaxidxs = lSortedGaps;

return;