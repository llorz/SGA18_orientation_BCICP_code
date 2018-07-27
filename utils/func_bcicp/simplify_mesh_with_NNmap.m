function [S_lowres, T_low2high, T_high2low] = simplify_mesh_with_NNmap(S,num_vtx,lowres_cache)
if nargin < 2, num_vtx = 1e3; end
if nargin < 3 
%% simplify the mesh 
S_lowres = struct();
surface = MESH.remesh(S.surface,struct('vertices',num_vtx));
S_lowres.surface.TRIV = surface.TRIV;
S_lowres.surface.VERT = surface.VERT;
S_lowres.surface.X = surface.VERT(:,1);
S_lowres.surface.Y = surface.VERT(:,2);
S_lowres.surface.Z = surface.VERT(:,3);
S_lowres.nf = surface.nf;
S_lowres.nv = surface.nv;
S_lowres.name = [S.name,'_',num2str(num_vtx/1e3),'k'];
% preprocess S_lowres
S_lowres = MESH.preprocess(S_lowres,'IfComputeLB',true,'IfFindNeigh',true,...
    'IfFindEdge',true,'IfComputeGeoDist',true);
else
    S_lowres = load([lowres_cache,S.name,'.mat']);
end
%% find the smooth nearest map between the original and simplified mesh
% nn-map: lowres-to-highres
T_low2high = knnsearch(S.surface.VERT, S_lowres.surface.VERT);
% nn-map: highres-to-lowres
T_high2low = knnsearch(S_lowres.surface.VERT, S.surface.VERT);

% if the NNmap is not smooth -> interpolate to remove outliers
v_id = get_discontinuous_vertex(T_low2high,S_lowres,S);
if ~isempty(v_id)
    T_low2high(v_id) = NaN;
    T_low2high = pMap_NNinterp(T_low2high,S_lowres);
end

v_id = get_discontinuous_vertex(T_high2low,S,S_lowres);
if ~isempty(v_id)
    T_high2low(v_id) = NaN;
    T_high2low = pMap_NNinterp(T_high2low,S);
end
end
%%
function [T12] = pMap_NNinterp(T12_in,S1)
s1_mapped = find(~isnan(T12_in));
s1_unmapped = find(isnan(T12_in));
[~, vid] = min(S1.Gamma(s1_unmapped, s1_mapped),[],2);
T12 = T12_in;
T12(s1_unmapped) = T12_in(s1_mapped(vid));
end


function [v_id, MAX_EDGE_LENGTH] = get_discontinuous_vertex(T21,S2,S1)
get_mat_entry = @(M,I,J) M(sub2ind(size(M),I,J));
dist1 = compute_mesh_dist_matrix(S1);
dist2 = compute_mesh_dist_matrix(S2);

d2 = get_mat_entry(dist2,S2.Elist(:,1),S2.Elist(:,2)); % edge length on S2
d1 = get_mat_entry(dist1,S1.Elist(:,1),S1.Elist(:,2)); % edge length on S1
MAX_EDGE_LENGTH = max(max(d2),max(d1));

d2_mapped = get_mat_entry(dist1, T21(S2.Elist(:,1)), T21(S2.Elist(:,2)));
% distortion = d2_mapped./d2; % edge distortion
e_id = d2_mapped > MAX_EDGE_LENGTH; % edges with large error
v_id = unique(S2.Elist(e_id,:)); % vertex with large error
v_id = reshape(v_id,[],1);
end