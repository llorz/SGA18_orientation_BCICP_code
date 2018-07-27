% 2018-05-10
% improve coverage + remove outliers + improve smoothness
function [T21_new, T12_new] = refine_pMap(T21_ini, T12_ini, S1,S2, num_iter)
if nargin < 5, num_iter = 4; end
compute_coverage = @(T) length(unique(T))/length(T);
T12_new = T12_ini;
T21_new = T21_ini;
n1 = S1.nv;
n2 = S2.nv;
%%
T12_new = pMapCoverage_smooth_many2one(T12_new,S1,S2);
T21_new = pMapCoverage_smooth_many2one(T21_new,S2,S1);
T21_new = pMapCoverage_uncoveredVtxHaveMultiPreImages_from_invMap(T21_new,S2,S1,T12_new);
T12_new = pMapCoverage_uncoveredVtxHaveMultiPreImages_from_invMap(T12_new,S1,S2,T21_new);
T12_new = pMapCoverage_NNwithMultiPreImg_to_uncoveredVtx_refined(T12_new,S1,S2,T21_new);
T21_new = pMapCoverage_NNwithMultiPreImg_to_uncoveredVtx_refined(T21_new,S2,S1,T12_new);
T12_new = get_smoothed_pMap(T12_new,S1,S2,1);
T21_new = get_smoothed_pMap(T21_new,S2,S1,1);

thres = 1e-4;
for iter = 1:num_iter
    cov1 = thres;
    while ((compute_coverage(T12_new) - cov1)/cov1 > thres)
        cov1 = compute_coverage(T12_new);
        T21_new = pMapCoverage_NNwithMultiPreImg_to_uncoveredVtx(T21_new,S2,S1,T12_new);
        T12_new = pMapCoverage_NNwithMultiPreImg_to_uncoveredVtx(T12_new,S1,S2,T21_new);
        
        T21_new = pMapCoverage_uncoveredVtxHaveMultiPreImages_from_invMap(T21_new,S2,S1,T12_new);
        T12_new = pMapCoverage_uncoveredVtxHaveMultiPreImages_from_invMap(T12_new,S1,S2,T21_new);
    end
    T12_new = get_smoothed_pMap(T12_new,S1,S2,1);
    T21_new = get_smoothed_pMap(T21_new,S2,S1,1);
    
    cov1 = thres;
    while ((compute_coverage(T12_new) - cov1)/cov1 > thres)
        cov1 = compute_coverage(T12_new);
        T12_new = pMapCoverage6(T12_new,S1,S2,T21_new);
        T21_new = pMapCoverage6(T21_new,S2,S1,T12_new);
        
        T21_new = pMapCoverage_uncoveredVtxHaveMultiPreImages_from_invMap(T21_new,S2,S1,T12_new);
        T12_new = pMapCoverage_uncoveredVtxHaveMultiPreImages_from_invMap(T12_new,S1,S2,T21_new);
    end
    T12_new = get_smoothed_pMap(T12_new,S1,S2,1);
    T21_new = get_smoothed_pMap(T21_new,S2,S1,1);
    
    if compute_coverage(T21_new) > 0.5 && compute_coverage(T12_new) > 0.5
        T21_new = fix_pMap_outliers(T21_new,S2,S1);
        T12_new = fix_pMap_outliers(T12_new,S1,S2);
    end
end

end