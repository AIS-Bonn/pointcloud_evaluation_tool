# Pointcloud Evaluation Tool

This tool computes the *Mean Map Entropy* and the *Mean Plane Variance* of a point cloud. 

Furthermore a pointcloud is generated that encodes the *Map Entropy* and *Plane Variance* in each point. 

Use the point cloud viewer of your choice to visualize these measures (e.g. pcl_viewer). 

Install
-------
```
mkdir build 
cd build
cmake .. 
make
```
Usage
-----
```
./mean_map_entropy path/to/pointcloud.pcd [-stepsize int] [-radius double] [-punishSolitaryPoints] [-minNeighbors int]
```
* **stepsize:** stepsize used iterating through the pointcloud (default: 1)
* **radius:** radius used for search for neighboring points (default: 0.3)
* **punishSolitaryPoints:** punishes points with bad values where the number of neighbors is under minNeighbors (default: false)
* **minNeighbors:** threshold (default: 15)

---

Please cite our paper if you use this tool. 
Razlaw, J., Droeschel, D., Holz, D., & Behnke, S. (2015, September). [Evaluation of registration methods for sparse 3D laser scans][1]. In Mobile Robots (ECMR), 2015 European Conference on (pp. 1-7). IEEE.
[1]: http://www.ais.uni-bonn.de/papers/ECMR_2015_Razlaw.pdf
