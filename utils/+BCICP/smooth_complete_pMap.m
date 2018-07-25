%% 2018-05-0
function [T12_new] = smooth_complete_pMap(T12_ini,S1,S2,num_iter)
% if length(unique(T12_ini))/length(T12_ini) < 0.6
%     T12_new = T12_ini;
%     fprintf('Skip the smoothing process.\n')
%     return;
% end
argmax = @(d) find(d == max(d), 1, 'first');
argmin = @(d) find(d == min(d), 1, 'first');

X1 = S1.surface.VERT; X2 = S2.surface.VERT;
N1 = S1.vtx_neigh; N2 = S2.vtx_neigh;
T12_in = T12_ini;
%%
for iter = 1:num_iter
%     prob_ind = get_discontinuous_vertex(T12_in,S1,S2);
    prob_ind = reshape(1:S1.nv,[],1); % smooth all the vtx
    if ~isempty(prob_ind)
        find_search_space = @(T12, ind) arrayfun(@(vid) unique(T12(cell2mat(N1(vid)))),...
            ind,'Uni',false);
        sp2 = find_search_space(T12_ini,prob_ind); % search space on S2
        sp2_neigh = cellfun(@(vids) unique(cell2mat(N2(vids))),sp2,'Uni',false);
        % remove the prob_ind from the search space
        
        
        W2 = S2.W;
        sp2_adj = cellfun(@(vids) full(W2(vids,vids)),sp2_neigh,'UniformOutput',false); % adjacency matrix
        sp2_conncomp = cellfun(@(adj) conncomp(adj), sp2_adj,'UniformOutput',false); % connected component
        sp2_compID = cellfun(@(ic) unique(ic),sp2_conncomp,'Uni',false);
        sp2_compSize = cellfun(@(comp,compID) histc(comp,compID),sp2_conncomp,sp2_compID,'UniformOutput',false);
        sp2_largestComp = cellfun(@(compID,compSize) compID(argmax(compSize)),sp2_compID,sp2_compSize,'Uni',false);
        sp2_final = cellfun(@(vids,ic,largestCompID) vids(ic == largestCompID),...
            sp2_neigh, sp2_conncomp, sp2_largestComp,'Uni',false);
        
        sp2_finalSize = cellfun(@(sp) length(sp),sp2_final);
        ind = prob_ind(sp2_finalSize > 2); % search space with size > 2
        T12 = T12_ini;
        d = X2(T12,:) - X1;
        ave_d = cell2mat(cellfun(@(neigh) mean(d(neigh,:),1),N1(ind),'UniformOutput',false));
        Y = X1(ind,:) + ave_d;
        
        tar_old = T12(ind);
        tar_new = cellfun(@(vid, sp) sp(knnsearch(X2(sp,:),Y(vid,:))),...
            num2cell(1:length(ind))',sp2_final(sp2_finalSize > 2)); % changed here
        
        T12_new = T12_in; T12_new(ind) = tar_new;
        T12_in = T12_new;
    else
        T12_new = T12_in;
        break
    end
end
end


