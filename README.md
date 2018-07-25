# Continuous and Orientation-preserving Correspondence via Functional Maps

This is a complete implementation for the paper "Continuous and Orientation-preserving Correspondence via Functinal Maps" by Jing Ren, Adrien Poulenard, Peter Wonka and Maks Ovsjanikov.


Main algorithms
------------------
- Compute a functional map with an orientation-preserving/reversing operator: `compute_fMap_regular_with_orientationOp(...)`
- Refine the point-wise maps: `bcicp_refine(...)`


Main parameters
------------------
- **k1(k2)**: the number of Eigen-basis used of mesh S1(S2)
- **numTimes**: the time-scale parameter to compute the WKS descriptors
- **skipSize**: the skip size of the computed WKS descriptors to save runtime
- **beta**: the weight for the orientation-preserving/reversing term
- **numIters**: the number of iterations for the BCICP refinement step

Note that, for all the tests in the paper, **k1 = k2 = 50**, **numIters = 10**. The rest three parameters were (sloppily) tuned on a ramdon shape pair in a dataset, then applied to all the rest shape pairs in this collection. Our choices of parameters of each tested datasets are (also specified in the example scripts):

- For **FAUST** dataset, we set the parameters: **numTimes = 100, skipSize = 20, beta = 1e-4**.
- For **TOSCA isometric** dataset, we set the parameters: **numTimes = 100, skipSize = 20, beta = 1**.
- For **TOSCA non-isometric** dataset, we set the parameters: **numTimes = 50, skipSize = 5, beta = 1**.

Examples
------------------
### Self-symmetric maps
The script `Example_selfSymm_Fig3.m` reproduces the Fig.3 of the paper: for a given shape, we use the orientation-reversing operator to compute its self-symmetric map. 

<img src="/figs/eg_selfSymm.png" width="600">

Note that here we used the **same** set of parameters to compute the self-symmetric maps. It would work much better if these parameters(especially **numTimes** and **beta**) are tuned per dataset. Alternatively, BCICP can be added to refine the self-symmetric maps as well.


### WKS (wave kernel signatures) initialization (on TOSCA dataset)
The script `Example_WKSini_Fig13.m` reproduces the Fig.13 of the paper and `Example_WKSini_Fig14.m` reproduces Fig.14: computing a map using +directOp/symmOp + BCICP.

<img src="/figs/WKSeg_Iso.png" height="150">                      <img src="/figs/WKSeg_nonIso.png" height="150">

### SEG (segmentation) initialization (on FAUST dataset)
The script `Example_SEGini.m` shows an example of non-symmetric segmentation initialization (compare the ICP and BCICP refinement). For the segmentation part, please refer to the paper 'Robust Structure-based Shape Correspondence' by Yanir Kleiman and Maks Ovsjanikov (and code: https://github.com/hexygen/structure-aware-correspondence).


Contact
------------------
Please let us know (jing.ren@kaust.edu.sa) if you have any question regarding the algorithm/paper (ฅωฅ*) or you find any bugs in the implementation (ÒωÓױ). 
