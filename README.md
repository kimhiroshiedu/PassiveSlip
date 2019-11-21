# PassiveSlip
Calculate passive slip and considering stress shadow

## subfunctions
* AutoCodingGmsh.m
    Auto code Gmsh script. It requires plate shape files (xyz) and edge line of model region.
* TriMeshGenerate.m
    Generate triB_*.txt files from m files exported by Gmsh.
* pick_interface_edge.m
    Test script for choosing the any coordinates on the map.
* Expand_Allsample.m
    Expand all compressed files, cha_*.mat and combine to a tcha.mat.
* maid2ia.m
    Convert result mat files to new version.
* out_passiveanalysis_vi.m
    Export variable files to make figures, and tables (for MCMC result files).
* out_passiveanalysis_v1.m
    Same as the above, but for forward simulation result.

## Plate data
* Meshes
    Mesh script and data files.
    * model_hirose
        * plate_hirose_pac.msh
        * plate_hirose_phs.msh
        * plate_hirose_phssagami.msh
        * plate_hirose_pac.m
        * plate_hirose_phs.m
        * plate_hirose_phssagami.m
        * tri_pac.txt
        * tri_phs.txt
        * tri_phssagami.txt
        * opt_phcont.txt
        Edge line of Nankai plate interface region.
        * opt_plate_phs.txt
        Depth data of Philippine Sea Plate interface of Nankai.
        * scale_0.25
            * plate_hirose_pac.msh
            * plate_hirose_phs.msh
            * plate_hirose_phssagami.msh
    * model_iwasaki
        * plate_iwasaki_pac.msh
        * plate_iwasaki_phs.msh
        * plate_iwasaki_phssagami.msh
        * plate_iwasaki_pac.m
        * plate_iwasaki_phs.m
        * plate_iwasaki_phssagami.m
        * tri_pac.txt
        * tri_phs.txt
        * tri_phssagami.txt
        * scale_0.15
            * plate_iwasaki_pac.msh
            * plate_iwasaki_phs.msh
            * plate_iwasaki_phssagami.msh
        * scale_0.25
            * plate_iwasaki_pac.msh
            * plate_iwasaki_phs.msh
            * plate_iwasaki_phssagami.msh
    * model_slab2 (now preparing...)
    * model_slab1 (now preparing...)