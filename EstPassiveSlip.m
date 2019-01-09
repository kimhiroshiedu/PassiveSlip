%% EstPassiveSlip
function EstPassiveSlip
% Coded by H. Kimura in 2019/1/9
warning('off','all')
input.parfile='./PARAMETER/parameter.txt';
input.optfile='./PARAMETER/opt_bound_par.txt';
% READ PARAMETER FOR MCMC Inversion 
[prm]=Read_Parameters(input);
% READ OBSERVATION FILE
[obs]=Read_obs(prm.file_obs);
% READ BLOCK BOUNDARY FILE in DIRECTORY
[BLK,obs]=READ_BLOCK_BOUND(prm.dir_block,obs);
% READ BLOCK INTERFACE BOUNDARY in DIRECTORY 
[BLK]=READ_BLOCK_INTERFACE(BLK,prm);
% SHOW BLOCK BOUNDARY MAP
SHOW_BLOCK_BOUND(BLK)
% READ FIX EULER POLES
[POL,prm]=READ_EULER_POLES(BLK,prm);
% READ RIGID BLOCK BOUNDARY
[BLK,prm]=READ_RIGID_BOUND(BLK,prm);
% READ INTERNAL DEFORMATION PARAMETER
[BLK,prm]=READ_INTERNAL_DEFORMATION(BLK,obs,prm);
% CALC. GREEN FUNCTION
[TRI,obs]=GREEN_TRI(BLK,obs);
% Combain to Green function
[D,G]=COMB_GREEN(BLK,obs,TRI);
% CALC. ABIC AND BLOCK MOTION
[BLK,obs]=CALC_AIC(BLK,obs,POL,prm);
% BLOCK MOTION BETWEEN TWO BLOCKS
[BLK,obs]=Est_Motion_BLOCKS(BLK,obs);
% MAKE FIGURES
MAKE_FIGS(BLK,obs);
% CAL Markov chain Monte Calro
[CHA]=MH_MCMC(D,G,BLK,prm,obs,POL);
% MAKE FIGURES
%MAKE_FIG(CHA,BLK,OBS,PRM);
% OUTPUT.DIR='./Result/';
WRITE_CHA(CHA,BLK,TRI,prm,obs,D,G)
%
end
%% Read Parameter file
function [prm]=Read_Parameters(input)
% Coded    by Hiroshi Kimura 2019/1/9 (ver 1.0)
%
fid=fopen(input.parfile,'r');
prm.home_d=pwd;
file_obs=fscanf(fid,'%s \n',[1,1]);
prm.file_obs=fullfile(prm.home_d,file_obs);
[~]=fgetl(fid);
dir_block=fscanf(fid,'%s \n',[1,1]);
prm.dir_block=fullfile(prm.home_d,dir_block);
[~]=fgetl(fid);
dir_block_interface=fscanf(fid,'%s \n',[1,1]);
prm.dir_block_interface=fullfile(prm.home_d,dir_block_interface);
[~]=fgetl(fid);
file_pole=fscanf(fid,'%s \n',[1,1]);
prm.file_pole=fullfile(prm.home_d,file_pole);
[~]=fgetl(fid);
file_rigb=fscanf(fid,'%s \n',[1,1]);
prm.file_rigb=fullfile(prm.home_d,file_rigb);
[~]=fgetl(fid);
file_internal=fscanf(fid,'%s \n',[1,1]);
prm.file_internal=fullfile(prm.home_d,file_internal);
[~]=fgetl(fid);
dir_result=fscanf(fid,'%s \n',[1,1]);
prm.dir_result=fullfile(prm.home_d,dir_result);
[~]=fgetl(fid);
%
prm.gpu=fscanf(fid,'%d \n',[1,1]);
[~]=fgetl(fid);
prm.itr=fscanf(fid,'%d \n',[1,1]);
[~]=fgetl(fid);
prm.thr=fscanf(fid,'%d \n',[1,1]);
[~]=fgetl(fid);
prm.cha=fscanf(fid,'%d \n',[1,1]);
[~]=fgetl(fid);
prm.kep=fscanf(fid,'%d \n',[1,1]);
[~]=fgetl(fid);
prm.rwd=fscanf(fid,'%f \n',[1,1]);
fclose(fid);
%====================================================
tmp=load(input.optfile);
prm.num=size(tmp,1);
prm.optB1=tmp(:,1);
prm.optB2=tmp(:,2);
prm.optINT=tmp(:,3);
%====================================================
fprintf('==================\nINPUT PARAMETERS\n==================\n') 
fprintf('Home directory           : %s \n',prm.home_d) 
fprintf('File observation data    : %s \n',prm.file_obs) 
fprintf('Directory block          : %s \n',prm.dir_block)
fprintf('Directory block interface: %s \n',prm.dir_block_interface) 
fprintf('File fixed Euler pole    : %s \n',prm.file_pole) 
fprintf('File Rigid boundary      : %s \n',prm.file_rigb) 
fprintf('File Internal deformation: %s \n',prm.file_internal) 
fprintf('Directory results        : %s \n',prm.dir_result) 
fprintf('GPU device (CPU:99)      : %i \n',prm.gpu) 
fprintf('Iteration (Max)          : %i \n',prm.itr) 
fprintf('Iteration (Threshold)    : %i \n',prm.thr) 
fprintf('Number of chain          : %i \n',prm.cha) 
fprintf('Number of keep           : %i \n',prm.kep) 
fprintf('Random walk distance     : %4.2f \n',prm.rwd) 
fprintf('==================\n') 
%====================================================
disp('PASS READ_PARAMETERS')
end
%% Read observation data
function obs=Read_obs(file_obs)
%-------------------
% INPUT format Observations:
% unit is mm/yr
% site_name lon lat EW_comp. NS_comp. UD_comp. ERR_EW ERR_NS ERR_UD 
%-------------------
fid_obs=fopen(file_obs,'r');
n=0;
while 1
  tline=fgetl(fid_obs);
  if ~ischar(tline); break; end
  str=strsplit(tline);
  n=n+1;
  obs(1).name(n) =cellstr(str(1));
  obs(1).alon(n) =str2double(cellstr(str(2))); %LON
  obs(1).alat(n) =str2double(cellstr(str(3))); %LAT
  obs(1).ahig(n) =str2double(cellstr(str(4))); %HIG
  obs(1).evec(n) =str2double(cellstr(str(5))); %E-W
  obs(1).nvec(n) =str2double(cellstr(str(6))); %N-S
  obs(1).hvec(n) =str2double(cellstr(str(7))); %U-D
  obs(1).eerr(n) =str2double(cellstr(str(8))); %E-W
  obs(1).nerr(n) =str2double(cellstr(str(9))); %N-S
  obs(1).herr(n) =str2double(cellstr(str(10))); %U-D
  obs(1).weight(n) =str2double(cellstr(str(11))); %Weight
  obs(1).axyz(n,:)=conv2ell(obs(1).alat(n),obs(1).alon(n),obs(1).ahig(n));
end
obs(1).nobs=n;
obs(1).ablk=zeros(obs(1).nobs,1);
obs(1).latmax=max(obs(1).alat);
obs(1).latmin=min(obs(1).alat);
obs(1).lonmax=max(obs(1).alon);
obs(1).lonmin=min(obs(1).alon);
fprintf('==================\n') 
fprintf('Number of GNSS site %4d \n',n)
fprintf('==================\n') 
end