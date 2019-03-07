%% Main part
function out_passiveanalysis_v1(folder)
result_dir = fullfile(pwd,folder);
fprintf('Now loading %s ...',fullfile(result_dir,'/obs.mat'))
load(fullfile(result_dir,'/obs.mat'));fprintf('load\n')
fprintf('Now loading %s ...',fullfile(result_dir,'/blk.mat'))
load(fullfile(result_dir,'/blk.mat'));fprintf('load\n')
fprintf('Now loading %s ...',fullfile(result_dir,'/cal.mat'))
load(fullfile(result_dir,'/cal.mat'));fprintf('load\n')
fprintf('Now loading %s ...',fullfile(result_dir,'/grn.mat'))
load(fullfile(result_dir,'/grn.mat'));fprintf('load\n')

ExportSlipDeficit(result_dir,cal,blk,d);
ExportAsperities(folder,blk);
ExportVectors(result_dir,obs,cal);
ExportStressStrain(result_dir,blk,cal);

fprintf('Files exported! \n');
end

%% Save slip deficit to txt files
function ExportSlipDeficit(folder,cal,blk,d)
mc = 1;
mr = 1;
dir_out=[folder,'/backslip'];
exid=exist(dir_out);
if exid~=7; mkdir(dir_out); end

for nb1=1:blk(1).nblock
  for nb2=nb1+1:blk(1).nblock
    nf=size(blk(1).bound(nb1,nb2).blon,1);
    if nf~=0
      fid     = fopen([dir_out,'/bs_',num2str(nb1),'_',num2str(nb2),'.txt'],'w');
      clon    = mean(blk(1).bound(nb1,nb2).blon,2);
      clat    = mean(blk(1).bound(nb1,nb2).blat,2);
      cdep    = mean(blk(1).bound(nb1,nb2).bdep,2);
      sdr     = sqrt( cal.slip(mc     :mc+  nf-1).^2 ...
                     + cal.slip(mc+  nf:mc+2*nf-1).^2 ...
                     + cal.slip(mc+2*nf:mc+3*nf-1).^2 );
      sdr_int = sdr.*d(1).id_lock(mc:mc+nf-1);
      flt_id  = mr:mr+nf-1;
      
      outdata = [flt_id' ...
          blk(1).bound(nb1,nb2).blon ...
          blk(1).bound(nb1,nb2).blat ...
          blk(1).bound(nb1,nb2).bdep ...
          clon clat cdep ...
          sdr sdr_int];
      
      fprintf(fid,'# %6s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %10s %10s \n',...
                      'tri_no','lon1','lon2','lon3','lat1','lat2','lat3','dep1','dep2','dep3',...
                      'c_lon','c_lat','c_dep','sdr[mm/yr]','sdr_int');
      fprintf(fid,'%8d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %10.4f %10.4f \n',outdata');
      fclose(fid);
      mr = mr+  nf;
      mc = mc+3*nf;
    end
  end
end

end

%% Export asperities
function ExportAsperities(folder,blk)
dir_out=[folder,'/backslip'];
exid=exist(dir_out);
if exid~=7; mkdir(dir_out); end
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    nf = size(blk(1).bound(nb1,nb2).blon,1);
    if nf ~= 0
      fid = fopen([dir_out,'/patchb_',num2str(nb1),'_',num2str(nb2),'.txt'],'w');
      for np = 1:size(blk(1).bound(nb1,nb2).patch,2)
        fprintf(fid,'> \n');
        patch = [blk(1).bound(nb1,nb2).patch(np).lon; ...
                 blk(1).bound(nb1,nb2).patch(np).lat];
        fprintf(fid,'%10.4f %10.4f \n',patch);
      end
      fclose(fid);
    end
  end
end

end

%% Export vectors of calculation, observation, residual vector
function ExportVectors(folder,obs,cal)
% 
% Calculate vector based on sampled parameter and green function
calvec.rig=cal.rig;
calvec.ela=cal.ela;
calvec.ine=cal.ine;
calvec.sum=cal.smp;
% Calculate residual vectors
obsv=[obs(1).evec;obs(1).nvec;obs(1).hvec];
calv=[cal.smp(1:3:end)';cal.smp(2:3:end)';cal.smp(3:3:end)'];
resvec.sum=obsv-calv;
% 
% Save to txt files
odir=[folder,'/vector'];
exid=exist(odir);
if exid~=7; mkdir(odir); end
obsvec=[obs(1).alon;obs(1).alat;obs(1).evec;obs(1).nvec;obs(1).hvec];
calvec.rig=[obs(1).alon;obs(1).alat;calvec.rig(1:3:end)';calvec.rig(2:3:end)';calvec.rig(3:3:end)'];
calvec.ela=[obs(1).alon;obs(1).alat;calvec.ela(1:3:end)';calvec.ela(2:3:end)';calvec.ela(3:3:end)'];
calvec.sum=[obs(1).alon;obs(1).alat;calvec.sum(1:3:end)';calvec.sum(2:3:end)';calvec.sum(3:3:end)'];
resvec.sum=[obs(1).alon;obs(1).alat;resvec.sum];
fid=fopen([odir,'/obs_vector.txt'],'w');
fprintf(fid,'%f %f %f %f %f\n',obsvec);
fclose(fid);
fid=fopen([odir,'/cal_vector.txt'],'w');
fprintf(fid,'%f %f %f %f %f\n',calvec.sum);
fclose(fid);
fid=fopen([odir,'/cal_vector_rig_site.txt'],'w');
fprintf(fid,'%f %f %f %f %f\n',calvec.rig);
fclose(fid);
fid=fopen([odir,'/cal_vector_ela_site.txt'],'w');
fprintf(fid,'%f %f %f %f %f\n',calvec.ela);
fclose(fid);
fid=fopen([odir,'/res_vector.txt'],'w');
fprintf(fid,'%f %f %f %f %f\n',resvec.sum);
fclose(fid);
% 
end

%% Export stress and strain
function ExportStressStrain(result_dir,blk,cal)
% Olivine
lambda = 74;
mu     = 82;
% 
mr = 1;
mc = 1;
odir = [result_dir,'/stress'];
exid = exist(odir);
if exid~=7; mkdir(odir); end
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    nf = size(blk(1).bound(nb1,nb2).blon,1);
    if nf ~= 0
      flt_id  = mr:mr+nf-1;
      clon    = mean(blk(1).bound(nb1,nb2).blon,2);
      clat    = mean(blk(1).bound(nb1,nb2).blat,2);
      cdep    = mean(blk(1).bound(nb1,nb2).bdep,2);
      e.xx = zeros(nf,1);
      e.xy = zeros(nf,1);
      e.yy = zeros(nf,1);
      e.xz = cal.strain(mc     :mc+  nf-1);
      e.yz = cal.strain(mc+  nf:mc+2*nf-1);
      e.zz = cal.strain(mc+2*nf:mc+3*nf-1);
      e_int.xx = zeros(nf,1);
      e_int.xy = zeros(nf,1);
      e_int.yy = zeros(nf,1);
      e_int.xz = cal.intstrain(mc     :mc+  nf-1);
      e_int.yz = cal.intstrain(mc+  nf:mc+2*nf-1);
      e_int.zz = cal.intstrain(mc+2*nf:mc+3*nf-1);
      [s]     = StrainToStress(e,lambda,mu);
      [s_int] = StrainToStress(e_int,lambda,mu);
      
      fid = fopen([odir,'/st_',num2str(nb1),'_',num2str(nb2),'.txt'],'w');
      outdata = [flt_id' ...
          blk(1).bound(nb1,nb2).blon ...
          blk(1).bound(nb1,nb2).blat ...
          blk(1).bound(nb1,nb2).bdep ...
          clon clat cdep ...
          e.xz e.yz e.zz ...
          e_int.xz e_int.yz e_int.zz ...
          s.xz s.yz s.zz ...
          s_int.xz s_int.yz s_int.zz];
      fprintf(fid,'# %6s %7s %7s %7s %7s %7s %7s %7s %7s %7s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s \n', ...
                  'tri_no','lon1','lon2','lon3','lat1','lat2','lat3','dep1','dep2','dep3',...
                  'c_lon','c_lat','c_dep', ...
                  'e_st','e_dp','e_ts','e_ini_st','e_ini_dp','e_ini_ts', ...
                  's_st','s_dp','s_ts','s_ini_st','s_ini_dp','s_ini_ts');
      fprintf(fid,'%8d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n',outdata');
      fclose(fid);

      mr = mr+  nf;
      mc = mc+3*nf;
    end
  end
end

end

%%
function Stress = StrainToStress(Strain, lambda, mu)
% StressToStrain.m
%
% Calculate stresses and invariants given a strain tensor and elastic
% moduli lambda and mu.
%
% Implements algorithms described in the journal article:
% Meade, B. J. Algorithms for calculating displacements, 
% strains, and stresses for triangular dislocation elements
% in a uniform elastic half space
% Computers and Geosciences, submitted, 2006.
%
% Use at your own risk and please let me know of any bugs/errors!
%
% Copyright (c) 2006 Brendan Meade
% 
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.

Stress.xx = 2.*mu.*Strain.xx + lambda.*(Strain.xx+Strain.yy+Strain.zz);
Stress.yy = 2.*mu.*Strain.yy + lambda.*(Strain.xx+Strain.yy+Strain.zz);
Stress.zz = 2.*mu.*Strain.zz + lambda.*(Strain.xx+Strain.yy+Strain.zz);
Stress.xy = 2.*mu.*Strain.xy;
Stress.xz = 2.*mu.*Strain.xz;
Stress.yz = 2.*mu.*Strain.yz;
Stress.I1 = Stress.xx + Stress.yy + Stress.zz;
Stress.I2 = -(Stress.xx.*Stress.yy + Stress.yy.*Stress.zz + Stress.xx.*Stress.zz) + Stress.xy.*Stress.xy + Stress.xz.*Stress.xz + Stress.yz.*Stress.yz;
Stress.I3 = Stress.xx.*Stress.yy.*Stress.zz + 2.*Stress.xy.*Stress.xz.*Stress.yz - (Stress.xx.*Stress.yz.*Stress.yz + Stress.yy.*Stress.xz.*Stress.xz + Stress.zz.*Stress.xy.*Stress.xy);

end
