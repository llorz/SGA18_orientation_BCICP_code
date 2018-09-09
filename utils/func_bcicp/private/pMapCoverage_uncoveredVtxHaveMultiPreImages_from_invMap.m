% Given T21: S2 -> S1 and T12, we want to use T12 to fix T21
% check if T21 of the uncovered vtx of T21 on S1 are mapped to some vtx
% with multiple pre-images from T21, if so, it's safe to remap those vtx.
function T21_new = pMapCoverage_uncoveredVtxHaveMultiPreImages_from_invMap(T21,S2,S1,T12_new)
T21_new = T21; n1 = S1.nv;
freq = get_pMap_range(T21_new,n1);
s1_ncid = reshape(setdiff(1:n1,T21_new),[],1); % vertices on S1 that are not covered
v_id = s1_ncid;
% to make the map bijective, we hope s_id can be mapped to v_id (given T12_new)
s_id = T12_new(v_id);
% if T21_new(s_id) has multiple preimages, we can remap those s_id
change_id = find(freq(T21_new(s_id)) > 1); % safe to remap these vertices
T21_new(s_id(change_id)) = v_id(change_id);
end