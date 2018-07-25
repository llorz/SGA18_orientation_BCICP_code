% 2018-06-11
clc; clear; close all;
addpath(genpath('utils/'))

k1 = 50;
k2 = 50;

mesh_dir = 'data/';
cache_dir = 'data/';

s1_name = 'tr_reg_000';
s2_name = 'tr_reg_004';


S1 = MESH.preprocess([mesh_dir, s1_name],'cacheDir',cache_dir,'numEigs',k1);
S2 = MESH.preprocess([mesh_dir, s2_name],'cacheDir',cache_dir,'numEigs',k2);
[S1.normals_vtx, S1.normals_face] = MESH.compute_vtx_and_face_normals(S1);
[S2.normals_vtx, S2.normals_face] = MESH.compute_vtx_and_face_normals(S2);
B1 = S1.evecs(:,1:k1); Ev1 = S1.evals(1:k1);
B2 = S2.evecs(:,1:k2); Ev2 = S2.evals(1:k2);
%% regular ini
t = 100;
fct1 = waveKernelSignature(B1, Ev1, S1.A, t);
fct2 = waveKernelSignature(B2, Ev2, S2.A, t);
fct1 = fct1(:,1:20:end);
fct2 = fct2(:,1:20:end);
%% ini + direct operator
para.beta = 1;
[C12_direct] = compute_fMap_regular_with_orientationOp(S1,S2,B1,B2,Ev1,Ev2,fct1,fct2,'direct',para);
T21_direct = annquery(C12_direct*B1', B2',1);

[C12_symm] = compute_fMap_regular_with_orientationOp(S1,S2,B1,B2,Ev1,Ev2,fct1,fct2,'symmetric',para);
T21_symm = annquery(C12_symm*B1', B2',1);
%%
plotOptions = {'OverlayAxis','y','cameraPos',[0,90]};
subplot(1,2,1);
MESH.PLOT.visualize_map_colors(S2,S1,T21_direct,plotOptions{:}); title('direct map')
subplot(1,2,2);
MESH.PLOT.visualize_map_colors(S2,S1,T21_symm,plotOptions{:}); title('symmetric map')