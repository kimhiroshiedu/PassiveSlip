# PassiveSlip
You can estimates Mechanical Coupling Distribution using MCMC, or do forward simulation by Est_PassiveSlip.m.  
Usage of main script, Est_PassiveSlip

     EST_PassiveSlip('option')

'option' is calculation method: 'fwd', forward simulation; 'mcmc', MCMC estimation  
If you enter no option (i.e. just type EST_PassiveSlip), MCMC estimation will start.
* Required files to run
    * Block line files in BLOCK directory.  
      File name requires block number, such as 01_hoge.txt, 02_fuga.txt  
      File format is Lon. and Lat. in 1st and 2nd columns, respectively. Note that the block line must be closed line (i.e. first and last rows must be equal).
    * (Option) Shape or trimesh files of Block boundaries.
    * Observation data file  
      File format is  

          Site_ID  Lon  Lat  Hig  vE  vN  vU  sigmaE  sigmaN  sigmaU  weight

      Units of each value are below  
      Lon, Lat: (deg)  
      Hig: (m)  
      vE, vN, vU, sigmaE, sigmaN, sigmaU: (mm/yr)
    * Parameter files
        * PARAMETER/Parameter.txt  
          Parameters to run the script.
        * PARAMETER/opt_bound_par.txt  
          Optimize the point interval of block boundary.
        * PARAMETER/interp_randwalkline.txt  
          Interpolate the udline in which the up- and down-dip limit (Zu and Zd) of mechanical coupling region.  
          File format is  

                Lon_Zd_1 Lat_Zd_1 Zd_1 Lon_Zu_1 Lat_Zu_1 Zu_1
                Lon_Zd_2 Lat_Zd_2 Zd_2 Lon_Zu_2 Lat_Zu_2 Zu_2
                                      •
                                      •
                                      •
                Lon_Zd_N Lat_Zd_N Zd_N Lon_Zu_N Lat_Zu_N Zu_N

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