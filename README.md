```
usage: g_aggregate.py [-h] -f TRAJ -s FILE [-o DIR] [-p STR] [-t NUMBER] [-n NUMBER] [--frames NUMBER]
                      [-m {residue,position}]

Script to identify TO clusters size and position along the MD trajectory.

optional arguments:
  -h, --help            show this help message and exit
  -f TRAJ, --traj TRAJ  Trajectory (all formats accepted by MDAnalysis)
  -s FILE, --topo FILE  Topology (all formats accepted by MDAnalysis)
  -o DIR, --outdir DIR  Output directory (default : .)
  -p STR, --prefix STR  Output prefix (default : Trajectory file name)
  -t NUMBER, --threshold NUMBER
                        Threshold for clustering (default : 13)
  -n NUMBER, --to-keep NUMBER
                        Number of largest clusters to keep for the results report (default : 2)
  --frames NUMBER       Number of frames to process (default : all)
  -m {residue,position}, --method {residue,position}
                        Method for clusters correspondance through frames (residue or position)
```
