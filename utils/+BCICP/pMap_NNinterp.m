% 2018-05-09
% T12_in: S1 -> S2 has NaN entries, i.e., no correspondences
% assign them the same target as their nearest neighbor
function [T12] = pMap_NNinterp(T12_in,S1)
s1_mapped = find(~isnan(T12_in));
s1_unmapped = find(isnan(T12_in));
[~, vid] = min(S1.Gamma(s1_unmapped, s1_mapped),[],2);
T12 = T12_in;
T12(s1_unmapped) = T12_in(s1_mapped(vid));
end