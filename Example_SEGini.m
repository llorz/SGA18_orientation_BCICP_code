clc; close all; clear;
addpath(genpath('utils/'))
%% parameters for FAUST
k1 = 50; k2 = 50; % #LB on source and target shape
numTimes = 100; % time-scale parameter to compute WKS descriptors
skipSize = 20; % skip size of the computed WKS descriptors - to save runtime
para.beta = 1e-1; % weight for the orientataion term
num_iters = 3; % num_iters for BCICP (in the paper we used 10)
plotOptions = {'IfShowCoverage',true,'cameraPos',[0,90],'OverlayAxis','y'};
%%
mesh_dir = 'data/FAUST/';
s1_name = 'tr_reg_002';
s2_name = 'tr_reg_085';
cache_dir = [mesh_dir,'cache/'];
if ~isdir(cache_dir), mkdir(cache_dir); end;
%% preprocess the mesh or load from the cache
S1 = MESH.preprocess([mesh_dir,s1_name], 'cacheDir',cache_dir,...
    'IfComputeLB',true, 'numEigs',k1,...
    'IfComputeGeoDist',true,...
    'IfComputeNormals',true);

S2 = MESH.preprocess([mesh_dir,s2_name], 'cacheDir',cache_dir,...
    'IfComputeLB',true, 'numEigs',k2,...
    'IfComputeGeoDist',true,...
    'IfComputeNormals',true);

%% compute the SEG descriptors
B1 = S1.evecs(:,1:k1); Ev1 = S1.evals(1:k1);
B2 = S2.evecs(:,1:k2); Ev2 = S2.evals(1:k2);

% compute the structure-based segmentation
[seg1,seg2] = MESH.SEG.compute_shape_pair_segmentation(S1,S2,'cacheDir',cache_dir,...
    'numComponentsRange',[8,7]);
% ignore the unmatched regions (with label 0)
unique_seg_id = setdiff(intersect(unique(seg1),unique(seg2)),0);
get_region = @(seg_id) cellfun(@(s_id) find(seg_id == s_id),num2cell(unique_seg_id),'UniformOutput',false);

fct1 = fMAP.compute_descriptors(S1,'Regions',get_region(seg1),'numSkip',skipSize);
fct2 = fMAP.compute_descriptors(S2,'Regions',get_region(seg2),'numSkip',skipSize);
%% use the regular fMap pipeline since the segmentation is non-symmetric
C12_ini = fMAP.compute_fMap_regular(S1,S2,B1,B2,Ev1,Ev2,fct1,fct2);
C21_ini = fMAP.compute_fMap_regular(S2,S1,B2,B1,Ev2,Ev1,fct2,fct1);
T12_ini = fMAP.fMap2pMap(B2,B1,C21_ini);
T21_ini = fMAP.fMap2pMap(B1,B2,C12_ini);
%% use ICP for refinement
C21 = fMAP.icp_refine(B2,B1,C21_ini,num_iters);
T12_icp = fMAP.fMap2pMap(B2,B1,C21);

C12 = fMAP.icp_refine(B1,B2,C12_ini,num_iters);
T21_icp = fMAP.fMap2pMap(B1,B2,C12);
%% use BCICP for refinement
[T21_new, T12_new] =  bcicp_refine(S1,S2,B1,B2,T21_ini, T12_ini, num_iters);
%% visualize the maps
% the black regions show the vertices are NOT covered by the computed map
figure(1);
subplot(1,3,1); cov = MESH.PLOT.visualize_map_colors(S1,S2,T12_ini,plotOptions{:}); title({'SEGini',['Coverage = ',num2str(cov,'%.2f'),'%']});
subplot(1,3,2); cov = MESH.PLOT.visualize_map_colors(S1,S2,T12_icp,plotOptions{:}); title({'SEGini + ICP',['Coverage = ',num2str(cov,'%.2f'),'%']});
subplot(1,3,3); cov = MESH.PLOT.visualize_map_colors(S1,S2,T12_new,plotOptions{:}); title({'SEGini + BCICP',['Coverage = ',num2str(cov,'%.2f'),'%']});

figure(2);
subplot(1,3,1); cov = MESH.PLOT.visualize_map_colors(S2,S1,T21_ini,plotOptions{:}); title({'SEGini',['Coverage = ',num2str(cov,'%.2f'),'%']});
subplot(1,3,2); cov = MESH.PLOT.visualize_map_colors(S2,S1,T21_icp,plotOptions{:}); title({'SEGini + ICP',['Coverage = ',num2str(cov,'%.2f'),'%']});
subplot(1,3,3); cov = MESH.PLOT.visualize_map_colors(S2,S1,T21_new,plotOptions{:}); title({'SEGini + BCICP',['Coverage = ',num2str(cov,'%.2f'),'%']});

