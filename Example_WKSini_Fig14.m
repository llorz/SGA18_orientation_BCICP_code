clc; close all; clear;
addpath(genpath('utils/'))
%% parameters for TOSCA-Isometric
k1 = 50; k2 = 50; % #LB on source and target shape
numTimes = 100; % time-scale parameter to compute WKS descriptors
skipSize = 20; % skip size of the computed WKS descriptors - to save runtime
para.beta = 1e-1; % weight for the orientataion term
num_iters = 3; % num_iters for BCICP (in the paper we used 10)
plotOptions = {'IfShowCoverage',false,'cameraPos',[-30,10]};
%%
mesh_dir = 'data/TOSCA_Isometric/';
s1_name = 'cat0';
s2_name = 'cat2';

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

%% compute the WKS descriptors
B1 = S1.evecs(:,1:k1); Ev1 = S1.evals(1:k1);
B2 = S2.evecs(:,1:k2); Ev2 = S2.evals(1:k2);

fct1 = waveKernelSignature(B1, Ev1, S1.A, numTimes);
fct2 = waveKernelSignature(B2, Ev2, S2.A, numTimes);
fct1 = fct1(:,1:skipSize:end);
fct2 = fct2(:,1:skipSize:end);
%% Step 01: using orientation term to compute the initial direct/symm maps
% with the direct operator
[C12_direct] = compute_fMap_regular_with_orientationOp(S1,S2,B1,B2,Ev1,Ev2,fct1,fct2,'direct',para);
[C21_direct] = compute_fMap_regular_with_orientationOp(S2,S1,B2,B1,Ev2,Ev1,fct2,fct1,'direct',para);
T12_direct = fMAP.fMap2pMap(B2,B1,C21_direct);
T21_direct = fMAP.fMap2pMap(B1,B2,C12_direct);
% with the symmetric operator
[C12_symm] = compute_fMap_regular_with_orientationOp(S1,S2,B1,B2,Ev1,Ev2,fct1,fct2,'symmetric',para);
[C21_symm] = compute_fMap_regular_with_orientationOp(S2,S1,B2,B1,Ev2,Ev1,fct2,fct1,'symmetric',para);
T12_symm = fMAP.fMap2pMap(B2,B1,C21_symm);
T21_symm = fMAP.fMap2pMap(B1,B2,C12_symm);
%% Step 02: add the BCICP refinement to improve the maps
tic
[T21_direct_new, T12_direct_new] =  bcicp_refine(S1,S2,B1,B2,T21_direct, T12_direct, num_iters);
[T21_symm_new, T12_symm_new] =  bcicp_refine(S1,S2,B1,B2,T21_symm, T12_symm, num_iters);
toc
%% Compare to ICP
C21 = fMAP.icp_refine(B2,B1,C21_direct,num_iters);
T12_icp = fMAP.fMap2pMap(B2,B1,C21);

C21_symm = fMAP.icp_refine(B2,B1,C21_symm,num_iters);
T12_symm_icp = fMAP.fMap2pMap(B2,B1,C21_symm);
%% visualize the maps
figure(1);
subplot(1,3,1); cov = MESH.PLOT.visualize_map_colors(S1,S2,T12_direct,plotOptions{:}); title({'+ directOp',['Coverage = ',num2str(cov,'%.2f'),'%']});
subplot(1,3,2); cov = MESH.PLOT.visualize_map_colors(S1,S2,T12_icp,plotOptions{:}); title({'+ directOp + ICP',['Coverage = ',num2str(cov,'%.2f'),'%']});
subplot(1,3,3); cov = MESH.PLOT.visualize_map_colors(S1,S2,T12_direct_new,plotOptions{:}); title({'+ directOp + BCICP',['Coverage = ',num2str(cov,'%.2f'),'%']});

figure(2);
subplot(1,3,1); cov = MESH.PLOT.visualize_map_colors(S1,S2,T12_symm,plotOptions{:}); title({'+ symmOp',['Coverage = ',num2str(cov,'%.2f'),'%']});
subplot(1,3,2); cov = MESH.PLOT.visualize_map_colors(S1,S2,T12_symm_icp,plotOptions{:}); title({'+ symmOp + ICP',['Coverage = ',num2str(cov,'%.2f'),'%']});
subplot(1,3,3); cov = MESH.PLOT.visualize_map_colors(S1,S2,T12_symm_new,plotOptions{:}); title({'+ symmOp + BCICP',['Coverage = ',num2str(cov,'%.2f'),'%']});
