% If you want to compute self-symmetric map,
% please use the function "compute_selfSymmetric_fMap"
%%
clc; close all; clear;
addpath(genpath('utils/'));
%% parameters
% These parameters are tested on the following dataset:
% SHREC19, FAUST, TOSCA human and animal shapes
% and it gives reasonable results: the initial map gives the correct
% orientation and after some BCICP refinement, we can get some accurate
% self-symmetric map.

% Actually, you can try to increase k (more basis functions) and numTimes (more descriptors)
% to get more accurate initial maps
% If you find the final map is not in a good quality (even after BCICP),
% you can try to increase para.d, i.e., add more weight to the orientation term.


k = 10; % the number of basis used
numTimes = 1; % the time-scale to compute WKS descriptors
% weights for different energy terms
para.a = 1; % Descriptors preservation
para.b = 1; % Commutativity with descriptors
para.c = 1; % Commutativity with Laplacian 
para.d = 1; % orientation-reversing term

% options to preprocess
mesh_options = {'IfComputeLB',true,'numEigs',k,... % compute k LB basis
    'IfComputeNormals',true,... % compute vtx normals for orientation term
    'IfComputeGeoDist',false};  % do not compute the geodesic distance matrix
%%
mesh_dir = 'data/selfSymm/'; 
s_name = 'baby.off';

S = MESH.preprocess([mesh_dir,s_name],mesh_options{:}); % preprocess the mesh
B = S.evecs(:,1:k); Ev = S.evals(1:k);
fct = waveKernelSignature(B, Ev, S.A, numTimes); % compute the WKS descriptors

[C_symm] = compute_selfSymmetric_fMap(S,S,B,B,Ev,Ev,fct,fct,para); 
C_symm = fMAP.icp_refine(B,B,C_symm,5);  % add 5 iteration of ICP
T_symm = fMAP.fMap2pMap(B,B,C_symm); % convert to pointwise map
MESH.PLOT.visualize_map_colors(S,S,T_symm,'IfShowCoverage',false,plotOptions{:});  % visualize the map
return
%% use BCICP to refine the self-symmetric map
S = MESH.preprocess(S, 'IfComputeLB',true,'numEigs',50,...
    'IfComputeGeoDist',true,...
    'IfFindEdge',true);
T_refined =  bcicp_refine(S, S, S.evecs, S.evecs,T_symm, T_symm, 5);
MESH.PLOT.visualize_map_colors(S,S,T_refined,'IfShowCoverage',false,plotOptions{:});  % visualize the map


