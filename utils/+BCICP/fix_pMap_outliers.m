function [T21_new,ic_id] = fix_pMap_outliers(T21,S2,S1)
thres = 0.08*S2.nv; % #vtx_in_one_component < thres: save to fix
W = S2.W; % adjacency matrix for mesh S2.
[v_id, MAX_EDGE_LENGTH] = get_discontinuous_vertex(T21,S2,S1);
dist2 = compute_mesh_dist_matrix(S2);

% break the connectivies of the problematic vertices
W(v_id,:) = 0;
W(:,v_id) = 0;

% find the connected components after removing the erroneous vertices
ic = conncomp(W); % component ID
% find the largest connected component
[ic_unique,ic_count,ic_id] = freqency_table(ic,false);


vid_err = cell2mat(ic_id(ic_count < thres));
vid_correct = ic_id{1}; % vtx in the largest connected component
[nn_dist, nn_id] = min(dist2(vid_err,vid_correct),[],2);
% if too far away from the nearest_neigh -> do not update
change_vid = find(nn_dist < MAX_EDGE_LENGTH);

tmp = ic_id{1};

T21_new = T21;
T21_new(vid_err(change_vid)) = T21(vid_correct(nn_id(change_vid)));
end
%
%% use the old version
% function [T12_new] = fix_pMap_outliers(T12,S1,S2,T12_ini)
% 
% T12_new = T12;
% if ~isfield(S1,'Gamma')
%     dist = calc_dist_matrix(S1);
% else
%     dist = S1.Gamma;
% end
% [vid_search,vid_err] = get_vid_with_continuous_correspondence(T12,S1,S2);
% [~,ic] = min(dist(vid_err,vid_search),[],2);
% 
% if length(vid_err) > 0.8*length(vid_search) % do not correct the outliers
%     if nargin > 3
%         T12_new = T12_ini;
%     else
%         T12_new = T12;
%     end
% else % it's safe? to correct the outliers
%     if nargin >3
%         T12_new(vid_err) = T12_ini(vid_err);
%     else
%         T12_new(vid_err) = T12(vid_search(ic));
%     end
% end
% end
% 
% function [v_id, v_id_prob] = get_vid_with_continuous_correspondence(T21,S2,S1)
% % T21: S2 -> S1
% get_pair_dist = @(X,id1,id2) sqrt(sum((X(id1,:) - X(id2,:)).^2,2));
% 
% if ~isfield(S2,'Elist'), S2.Elist = get_mesh_edge_list(S2); end
% if ~isfield(S1,'Elist'), S1.Elist = get_mesh_edge_list(S1); end
% 
% X1 = S1.surface.VERT;
% X2 = S2.surface.VERT;
% d2 = get_pair_dist(X2,S2.Elist(:,1),S2.Elist(:,2));
% d1 = get_pair_dist(X1,S1.Elist(:,1),S1.Elist(:,2));
% d2_mapped = get_pair_dist(X1,T21(S2.Elist(:,1)),T21(S2.Elist(:,2)));
% 
% W = S2.W; % adjacency matrix
% e_id = d2_mapped > max(max(d2),max(d1)); % edges with large error
% v_id = unique(S2.Elist(e_id,:)); % vertex with large error
% 
% % break the connectivies of the problematic vertices
% W(v_id,:) = 0;
% W(:,v_id) = 0;
% 
% % find the connected components after removing the erroneous vertices
% ic = conncomp(W); % component ID
% 
% % find the largest connected components
% ic_unique = unique(ic);
% count = histc(ic,unique(ic));
% [~, ic_id] = max(count);
% 
% v_id = find(ic == ic_unique(ic_id));   % vtx with acceptable correspondences
% v_id_prob = find(ic ~= ic_unique(ic_id)); % vtx with bad correspondence
% v_id = reshape(v_id,[],1);
% v_id_prob = reshape(v_id_prob,[],1);
% end
% 
% % get the edge list and the corresponding edge weight
% function [Elist,K] = get_mesh_edge_list(S)
% T = S.surface.TRIV;
% 
% I = [T(:,1);T(:,2);T(:,3)];
% J = [T(:,2);T(:,3);T(:,1)];
% 
% Elisto = [I,J];
% sElist = sort(Elisto,2);
% [Elist,~] = unique(sElist, 'rows');
% 
% if nargout > 1  % get edge weight
%     n = size(S.A,1); % num of vtx
%     ind = sub2ind([n,n],Elist(:,1),Elist(:,2));
%     if isfield(S,'W')
%         W = S.W;
%         K = -W(ind);  % true edge weights
%     else
%         W = cotLaplacian([S.surface.X, S.surface.Y, S.surface.Z], S.surface.TRIV);
%         K = -W(ind);
%     end
% end
% end
% 
% 
