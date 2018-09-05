# Continuous and Orientation-preserving Correspondence via Functional Maps

This is a complete implementation for the paper "Continuous and Orientation-preserving Correspondence via Functinal Maps" by Jing Ren, Adrien Poulenard, Peter Wonka and Maks Ovsjanikov.


Main algorithms
------------------
- Compute a functional map with an orientation-preserving/reversing operator: 
```
C12 = compute_fMap_regular_with_orientationOp(S1,S2,B1,B2,Ev1,Ev2,fct1,fct2,type)

% Input:
%   S1: the source mesh with the new basis B1, and the corresponding eigenvalues Ev1 (k1 Eigen-functions)
%   S2: the target mesh with the new basis B2, and the corresponding eigenvalues Ev2 (k2 Eigen-functions)
%   fct1: the descriptors of shape S1
%   fct2: the descriptors of shape S2
%   type: 'direct' or 'symmetric' (call the orientation preserving/reversing operator)
% Output:
%   C12: a functional map from S1 -> S2 (k2-by-k1 matrix)
```
- Refine the point-wise maps: 
```
[T21, T12] = bcicp_refine(S1,S2,B1,B2,T21_ini,T12_ini,num_iter)

% Input: 
%   S1/S2 with corresponding Eigen-functions B1/B2
%   Initial point-wise maps from both directions: T12_ini: S1 -> S2 and T21_ini: S2 -> S1
%   num_iter: number of iterations to run BCICP refinement
% Output:
%   T12, T21: the refined point-wise maps with better accuracy, smoothness, bijectivity and coverage
```

Main parameters
------------------
- **numTimes**: the time-scale parameter to compute the WKS descriptors
- **skipSize**: the skip size of the computed WKS descriptors to save runtime
- **k1(k2)**: the number of Eigen-basis used of mesh S1(S2)
- **beta**: the weight for the orientation-preserving/reversing term
- **numIters**: the number of iterations for the BCICP refinement step

Note that, for all the tests in the paper, **k1 = k2 = 50**, **numIters = 10** and **beta = 0.1**. The rest two parameters for **computing the WKS descriptors** were (sloppily) tuned on a ramdon shape pair in a dataset, then applied to all the rest shape pairs in this collection. Our choices of parameters of each tested datasets are (also specified in the example scripts):

- For **FAUST** dataset, we set the parameters: **numTimes = 100, skipSize = 10**.
- For **TOSCA isometric** dataset, we set the parameters: **numTimes = 100, skipSize = 20**.
- For **TOSCA non-isometric** dataset, we set the parameters: **numTimes = 50, skipSize = 10**.

Examples
------------------
### Self-symmetric maps
The script `Example_selfSymm_Fig3.m` reproduces the Fig.3 of the paper: for a given shape, we use the orientation-reversing operator to compute its self-symmetric map. 

<img align="center"  src="/figs/eg_selfSymm.png" width="600">

Note that here we used the **same** set of parameters to compute the self-symmetric maps. It would work much better if these parameters (especially **numTimes** and **beta**) are tuned per dataset. Moreover, BCICP can be added to refine the self-symmetric maps as well.


### WKS (wave kernel signatures) initialization (on TOSCA dataset)
The script `Example_WKSini_Fig14.m` and `Example_WKSini_Fig15.m` reproduces Fig.14-15 of the paper: computing a map using +directOp/symmOp + BCICP.

Note that for TOSCA non-isometric dataset, the isometry assupmtion fails. The orientation-preserving/reversing operator still works to some extent but much less effective than the cases of isometric shape pair. In this case, more Eigen-basis are needed to compute the WKS desriptors and more BCICP iterations are needed to refine the maps. 

<img src="/figs/WKSeg_Iso.png" height="150"> &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  <img src="/figs/WKSeg_nonIso.png" height="150">

### SEG (segmentation) initialization (on FAUST dataset)
The script `Example_SEGini.m` shows an example of non-symmetric segmentation initialization (compare the ICP and BCICP refinement). For the segmentation part, please refer to the paper "Robust Structure-based Shape Correspondence" by Yanir Kleiman and Maks Ovsjanikov (and code: https://github.com/hexygen/structure-aware-correspondence).


Contact
------------------
Please let us know (jing.ren@kaust.edu.sa) if you have any question regarding the algorithms/paper (ฅωฅ*) or you find any bugs in the implementation (ÒωÓױ). 
