This module is to recognize and characterize TO molecule clusters in molecular dynamic simulations. 
It uses MDAnalysis library. 

## Installation

Tested with python 3.8. Feel free to give feedback for other python versions. 

Required packages : 
* MDAnalysis 
* Matplotlib v3.2.2

We advice to use conda or venv environment. 
``` 
pip install -r requirements.txt 
```

You need to add module to your PYTHONPATH (pip package will come) :
```
export PYTHONPATH=$PYTHONPATH:<this_directory>/src/
```

## Usage

```
usage: g_aggregate.py [-h] -f FILE -s FILE [-o DIR] [-p STR] [--c-mol STR] [-t NUMBER] [-n NUMBER]
                      [-m {residue,position}] [-c NUMBER] [--frames NUMBER] [--m-mol STR] [--slice-size FLOAT]
                      [--slice-axis STR] [--process INT] [--membrane-axis {x,y,z}]

Script to identify TO clusters size and position along the MD trajectory.

optional arguments:
  -h, --help            show this help message and exit

General options:
  -f FILE, --traj FILE  Trajectory (all formats accepted by MDAnalysis)
  -s FILE, --topo FILE  Topology (all formats accepted by MDAnalysis)
  -o DIR, --outdir DIR  Output directory (default : .)
  -p STR, --prefix STR  Output prefix (default : Trajectory file name)
  --frames NUMBER       Number of frames to process (default : all)
  --process INT         Number of process to run for membrane computation (default : available cpus)

Membrane options:
  --m-mol STR           Membrane molecule (default : DOPC)
  --slice-size FLOAT    Membrane slice size (default : 10)
  --slice-axis STR      Axis to slice membrane (default : x)
  --membrane-axis {x,y,z}
                        Axis to determine membrane highest and lowest points (default : z)

Clustering options:
  --c-mol STR           Clusters molecule (default : TO)
  -t NUMBER, --threshold NUMBER
                        Threshold for clustering (default : 13)
  -n NUMBER, --to-keep NUMBER
                        Number of largest clusters to keep (default : 2)
  -m {residue,position}, --method {residue,position}
                        Method for clusters correspondance through frames (residue or position) (default :
                        residue)
  -c NUMBER, --nb-corr NUMBER
                        Number of largest clusters to take into account for clusters correspondance at time t
                        (default : all)

```

## Clusters correspondance
In dynamic trajectory, we have to identify at which clusters at time t-1 corresponds clusters at time t.
Clusters are initialized at frame 0, cluster 0 is assigned to largest cluster, cluster 1 to 2nd largest cluster and so on.
Then, there's two methods to identify clusters in next frames : 
- By residues : with this method, the cluster at time t-1 is assigned to the cluster with "highest residue correspondance" at time t. Highest residue correspondance is computed this way between two clusters : number of identical residues over size of the largest cluster. 
- By position : with this method, the cluster at time t-1 is assigned to the closest cluster at time t, based on distance between center of masses. This method is longer, because all center of masses have to be computed. 

You can configure the methods by changing the number of clusters to take into account for clusters correspondance (-c/--nb-corr option). It means at time t, only the choosen number largest clusters will be compared to kept t-1 clusters (-n/-to-keep number). For example, if you choose all, all the clusters at time t will be compared with kept clusters at time t-1, and the closest will be assigned to each. Keep all is long for position method because all centers of mass have to be computed and the results will be biaised by membrane movement. 

## Outputs 

* \<prefix>_size_clusters.svf : Size of largest clusters through time
* \<prefix>_absolute_com.svg : Position of largest clusters raw center of masses through all frames. 
* \<prefix>_relative\_<z|x|y>.svg : x|y|z mean coordinates of clusters and highest and lowest point of membrane.



