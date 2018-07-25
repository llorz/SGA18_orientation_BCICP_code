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

Note that, for all the tests in the paper, **k1 = k2 = 50**, **numIters = 10**. The rest three parameters were (sloppily) tuned on a ramdon shape pair in a dataset, then applied to all the rest shape pairs in this collection. Our choices of parameters of each tested datasets are specified in the example scripts.


Examples
------------------
### Self-symmetric maps
The script `Example_selfSymm_Fig3` reproduces the Fig.3 of the paper: for a given shape, we use the orientation-reversing operator to compute its self-symmetric map. 

![Self-symmetric maps](/figs/eg_selfSymm.png = 800x)

Note that here we used the **same** set of parameters to compute the self-symmetric maps. It would work much better if these parameters(especially **numTimes** and **beta**) are tuned per dataset. Alternatively, BCICP can be added to refine the self-symmetric maps as well.


### WKS (wave kernel signatures) initialization (on TOSCA (non-)isometric dataset)
![TOSCA isometric pair](/figs/WKSeg_Iso.png)


### SEG (segmentation) initialization (on FAUST dataset)
![TOSCA non-isometric pair](/figs/WKSeg_nonIso.png)


Contact
------------------
