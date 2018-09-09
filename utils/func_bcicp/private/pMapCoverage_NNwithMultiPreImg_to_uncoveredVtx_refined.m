% for an uncovered vertex v on shape S1 wiht neighbors neigh (on shape S1)
% if some of the neighbors have multiple preimages, search among those
% preimages, and find the nearest neighbor w of T2(v),then it's safe to set
% T1(w) = v without hurting the coverage rate.
% a refined version of tmp_improve_fMap_bijectivity2.m
function T21_new =  pMapCoverage_NNwithMultiPreImg_to_uncoveredVtx_refined(T21,S2,S1,T12_new)
X2 = S2.surface.VERT;
get_mat_entry = @(M,I,J) M(sub2ind(size(M),I,J));
dist2 = compute_mesh_dist_matrix(S2);

N1 = S1.vtx_neigh;
n1 = S1.nv;

T21_new = T21;
freq = get_pMap_range(T21,n1);
nc_id_all = setdiff(1:n1,unique(T21_new))'; % all the vertices that are not covered by T1
s_id_all = cellfun(@(neigh) find(ismember(T21_new,neigh(freq(neigh) > 1))),...
    N1(nc_id_all),'UniformOutput',false);
change_id = find(cellfun(@(x) length(x),s_id_all) > 0);
nc_id_all_new = nc_id_all(change_id);
s_id_all_new = s_id_all(change_id);

new_preimage = cellfun(@(s_id,nc_id)s_id(knnsearch(X2(s_id,:),X2(T12_new(nc_id),:))),...
    s_id_all_new,num2cell(nc_id_all_new));

v_id = find(get_mat_entry(dist2,T12_new(nc_id_all_new),new_preimage) < 1e-3);
T21_new(new_preimage(v_id)) = nc_id_all_new(v_id);
end