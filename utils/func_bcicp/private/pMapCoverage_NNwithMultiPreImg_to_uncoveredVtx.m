% T21: S2 -> S1 (X2 -> X1)
% T12: S1 -> S2 (X1 -> X2)
% for an uncovered vertex v on shape S1, find the nearest neighbor w of T12(v)
% such that the T21(w) has multiple pre-images, then it's safe to set
% T1(w) = v without hurting the coverage rate.
function T21_new = pMapCoverage_NNwithMultiPreImg_to_uncoveredVtx(T21,S2,S1,T12)
get_pair_dist = @(X,id1,id2) sum((X(id1,:)- X(id2,:)).^2,2);
get_mat_entry = @(M,I,J) M(sub2ind(size(M),I,J));
dist1 = compute_mesh_dist_matrix(S1);
dist2 = compute_mesh_dist_matrix(S2);


X1 = S1.surface.VERT;
X2 = S2.surface.VERT;
n1 = S1.nv;
T21_new = T21;

% vtx ID: not covered by the pMap T1; on shape S1
nc_id = setdiff(1:n1,unique(T21));
% the range of pMap T1; on shape S1
freq = get_pMap_range(T21,n1);

% vtx ID: the vertices on S2 that their images(the vertices on S1) have
% multiple pre-images, i.e., it's "safe" to change their images
search_v = find(freq(T21_new) > 1); % on S2
src_id = T12(nc_id); % on S2
% nn_id = search_v(knnsearch(X2(search_v,:),X2(src_id,:))); % on S2
[~,ic] = min(dist2(src_id, search_v),[],2);
nn_id = search_v(ic);
if ~isempty(nn_id)
    % replace_id = find(abs(T1(T2(nn_id)) - nn_id) > abs(T1(nc_id) - nn_id));
    err_old = get_mat_entry(dist2,T12(T21(nn_id)),nn_id);
    err_new = get_mat_entry(dist2,T12(nc_id),nn_id);
    % err_old = get_pair_dist(X2,T12(T21(nn_id)),nn_id);
    % err_new = get_pair_dist(X2,T12(nc_id),nn_id);
    
    replace_id = find((err_old - err_new)./err_old > -2e-2);
    T21_new(nn_id(replace_id)) = nc_id(replace_id);
end
%% better but much slower: search_v is updated in each iteration
% for i = 1:length(nc_id)
%     src_id = T2(nc_id(i));
%     search_v = find(freq(T1_new) > 1);
%     nn_id = search_v(knnsearch(X2(search_v,:),X2(src_id,:)));
%     if (abs(T1(T2(nn_id)) - nn_id) > abs(T1(nc_id(i)) - nn_id))
%         T1_new(nn_id) = nc_id(i);
%     end
% end
end