%% Main part
function out_passiveanalysis_v1(folder)
result_dir = fullfile(pwd,folder);
fprintf('Now loading %s ...',fullfile(result_dir,'/prm.mat'))
load(fullfile(result_dir,'/prm.mat'));fprintf('load\n')
fprintf('Now loading %s ...',fullfile(result_dir,'/obs.mat'))
load(fullfile(result_dir,'/obs.mat'));fprintf('load\n')
fprintf('Now loading %s ...',fullfile(result_dir,'/blk.mat'))
load(fullfile(result_dir,'/blk.mat'));fprintf('load\n')
fprintf('Now loading %s ...',fullfile(result_dir,'/tri.mat'))
load(fullfile(result_dir,'/tri.mat'));fprintf('load\n')
fprintf('Now loading %s ...',fullfile(result_dir,'/grn.mat'))
load(fullfile(result_dir,'/grn.mat'));fprintf('load\n')
G(1).tb=full(G(1).tb);

ExportSlipDeficit(DIR,TCHA,BLK,SDR);
ExportVectors(DIR,BLK,TCHA,G,D,GRD,TRIg,OBS);
end

%% Save slip deficit to txt files
function ExportSlipDeficit(DIR,TCHA,BLK,SDR)
NN=1;
folder=[DIR,'/coupling'];
exid=exist(folder);
if exid~=7; mkdir(folder); end
FIDstdinfo=fopen([folder,'/Std_info.txt'],'w');
fprintf(FIDstdinfo,'NB1 NB2 STDmax STDmin\n');
for NB1=1:BLK(1).NBlock
  for NB2=NB1+1:BLK(1).NBlock
    NF=size(BLK(1).BOUND(NB1,NB2).blon,1);
    if NF~=0
      FIDmain = fopen([folder,'/C_',num2str(NB1),'_',num2str(NB2),'.txt'],'w');
      FLTNUM = NN:NN+NF-1;
      AVECP = TCHA.AVEFLT(FLTNUM,:);
      MEDCP = TCHA.MEDFLT(FLTNUM,:);
      SDRs  =    SDR.flax(FLTNUM,:);
      STD   = TCHA.STDFLT(FLTNUM,:);
      clon = mean(BLK(1).BOUND(NB1,NB2).blon,2);
      clat = mean(BLK(1).BOUND(NB1,NB2).blat,2);
      cdep = mean(BLK(1).BOUND(NB1,NB2).bdep,2);
      outdata = [FLTNUM' ...
          BLK(1).BOUND(NB1,NB2).blon ...
          BLK(1).BOUND(NB1,NB2).blat ...
          BLK(1).BOUND(NB1,NB2).bdep ...
          clon clat cdep ...
          AVECP MEDCP SDRs STD];
      fprintf(FIDmain,'# %6s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %10s %10s %10s %10s\n',...
                      'Tri_No','Lon1','Lon2','Lon3','Lat1','Lat2','Lat3','Dep1','Dep2','Dep3',...
                      'C_Lon','C_Lat','C_Dep','Mean_Cp','Median_Cp','SDR[mm/yr]','sigma');
      fprintf(FIDmain,'%8d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %10.4f %10.4f %10.4f %10.4f\n',outdata');
      fprintf(FIDstdinfo,'%d %d %f %f\n',NB1,NB2,min(STD),max(STD));
      fclose(FIDmain);
      sdr(1).bound(NB1,NB2).str=SDR.str(FLTNUM,:);
      sdr(1).bound(NB1,NB2).tns=SDR.tns(FLTNUM,:);
      sdr(1).bound(NB1,NB2).dip=SDR.dip(FLTNUM,:);
      NN=NN+NF;
    end
  end
end
fclose(FIDstdinfo);
sdrfile=[folder,'/sdr'];
save(sdrfile,'sdr','-v7.3')
end

%% Export vectors of calculation, observation, residual vector
function ExportVectors(DIR,BLK,TCHA,G,D,GRD,TRIg,OBS)
% 
calvec=calc_sampling_vector(OBS,BLK,TCHA,D,G);
resvec=calc_residual_vector(BLK,OBS,calvec);
[grdvec,GRD]=calc_vector_atmesh(BLK,TCHA,D,G,GRD,TRIg);
% 
WRITE_VECTOR(DIR,OBS,BLK,calvec,resvec,grdvec,GRD);
end
