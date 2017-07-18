PoseTrack: Joint Multi-Person Pose Estimation and Tracking
========================
 
# Introduction
The repo provides the source code for our CVPR'17 [paper](https://arxiv.org/abs/1611.07727) for joint mult-person pose estimation and tracking. 

<p align="left">
<img src="http://pages.iai.uni-bonn.de/iqbal_umar/PoseTrack/data/PoseTrack.gif", width="1000">
</p>

**Umar Iqbal, Anton Milan, Juergen Gall**  
PoseTrack: Joint Multi-Person Pose Estimation and Tracking  
IEEE Conference on Computer Vision and Pattern Recongnition (CVPR) 2017.  
[Project Page](http://pages.iai.uni-bonn.de/iqbal_umar/PoseTrack/)  

The code is tested on Ubuntu 16.04 (64bit) with MATLAB (2016a).  

## Installation

##### Dependencies
- C++ 11
- CUDA >=7.5
- MATLAB
- HDF5 1.8
- CMake
- [Caffe](http://caffe.berkeleyvision.org/installation.html)
- [Gurobi 6.0.x](https://user.gurobi.com/download/gurobi-optimizer)

##### Installation Instructions
1. Clone repository	
   ```
   $ git clone https://github.com/iqbalu/PoseTrack-CVPR2017.git --recursive
   ```

2. Build and download dependencies of [DeeperCut](https://github.com/eldar/deepcut)
..1. Build Caffe and its MATLAB interface after configuring Makefile.config
   ```
   $ cd deepcut/external/caffe
   $ make -j 12 all matcaffe
   ```

..2. Build `liblinear`, specify the path to the MATLAB installation	
   ```
   $ cd deepcut/external/liblinear-1.94/matlab
   $ CC=gcc CXX=g++ MATLABDIR=PATH_TO_MATLAB make
   ```

..3. Download models
   ```
   $ cd deepcut/data
   $ ./download_models.sh
   ```

3. Build PoseTrack solver	
   ```
   $ cd external/solver
   $ cmake . -DGUROBI_ROOT_DIR=/path/to/gurobi605/linux64 -DGUROBI_VERSION=60
   $ make solver-callback
   ```

4. Obtain Gurobi license from http://www.gurobi.com/downloads/licenses/license-center
   and place the license at desired location and modify the p.gurobi_license_file in posetrack/lib/pt_exp_params.m to point to the license file location.  

5. Download PoseTrack models
```
$ cd data
$ ./download_models.sh
```

## Run Demo	
```
% in MATLAB
>> startup
>> demo_posetrack
```

## Citing
```
@inproceedings{iqbal2016PoseTrack,
	author = {Umar Iqbal and Anton Milan and Juergen Gall},
	title = {PoseTrack: Joint Multi-Person Pose Estimation and Tracking},
	booktitle = {IEEE Conference on Computer Vision and Pattern Recognition (CVPR)},
	year = {2017},
	url = {https://arxiv.org/abs/1611.07727}
}
```

## Acknowledgements  
The authors are thankful to Chau Minh Triet, Andreas Doering, and Zain Umer Javaid for the help with annotating the dataset. The work has been financially supported by the DFG project GA 1927/5-1 (DFG Research Unit FOR 2535 Anticipating Human Behavior) and the ERC Starting Grant ARCA (677650).  
Thanks to Eldar Insafutdinov for releasing the code for [DeeperCut](https://github.com/eldar/deepcut).






