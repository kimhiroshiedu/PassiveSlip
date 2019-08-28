%% Main part
function out_passiveanalysis_vi(folder)
% for inversion results
savedir = fullfile(pwd,folder);
fprintf('Now loading %s ...',fullfile(savedir,'/prm.mat'))
load(fullfile(savedir,'/prm.mat'));fprintf('load\n')
fprintf('Now loading %s ...',fullfile(savedir,'/obs.mat'))
load(fullfile(savedir,'/obs.mat'));fprintf('load\n')
fprintf('Now loading %s ...',fullfile(savedir,'/blk.mat'))
load(fullfile(savedir,'/blk.mat'));fprintf('load\n')
fprintf('Now loading %s ...',fullfile(savedir,'/grn.mat'))
load(fullfile(savedir,'/grn.mat'));fprintf('load\n')
fprintf('Now loading %s ...',fullfile(savedir,'/tcha.mat'))
load(fullfile(savedir,'/tcha.mat'));fprintf('load\n')
G(1).tb_kin = full(G(1).tb_kin);
G(1).tb_mec = full(G(1).tb_mec);

[bslip,bslipl,vec] = CalcOptimumValue(prm,tcha,G,d);
SavePoles(savedir,blk,tcha);
SaveBackslip(savedir,blk,tcha);
SaveVectors(savedir,tcha,obs,vec);
SaveInternalStrain(savedir,tcha);
ExportSlipDeficit(savedir,cal,blk,d);
ExportAsperities(folder,blk);
ExportVectors(savedir,obs,cal);
ExportStressStrain(savedir,blk,d,G,cal);

fprintf('Files exported. \n');
end

%% Save internal strains
function SaveInternalStrain(folder,tcha)

end

%% Save obs-, cal-, and res-vectors
function SaveVectors(folder,tcha,obs,vec)
% save vobs
% save vcal
% save vres
end

%% Save coupling and backslip
function SaveBackslip(folder,blk,tcha)

end

%% Save Euler vectors
function SavePoles(folder,blk,tcha)
fid = fopen(fullfile(folder,'/est_euler_pole.txt'),'wt');
fprintf(fid,'# %6s %8s %7s %8s %11s %9s %9s %9s %9s %9s %9s\n',...
    'Blk_no','Lat(deg)','Lon(deg)','Ang(deg/my)',...
    'sigxx','sigxy','sigxz','sigyy','sigyz','sigzz');
fprintf('# %6s %8s %7s %8s %11s %9s %9s %9s %9s %9s %9s\n',...
    'Blk_no','Lat(deg)','Lon(deg)','Ang(deg/my)',...
    'sigxx','sigxy','sigxz','sigyy','sigyz','sigzz');
for bk = 1:blk(1).nblock
  [latp,lonp,ang] = xyzp2lla(tcha.avepol(3.*bk-2,:),tcha.avepol(3.*bk-1,:),tcha.avepol(3.*bk,:));
  [a,b,c,d,e,f] = out_cov(tcha.covpol(3*bk-2:3*bk,3*bk-2:3*bk));
  fprintf('%8i %7.2f %8.2f %11.2e %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f \n',...
    bk,mean(latp),mean(lonp),mean(ang),a,b,c,d,e,f);
  fprintf(fid,'%8i %7.2f %8.2f %11.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f \n',...
    bk,mean(latp),mean(lonp),mean(ang),a,b,c,d,e,f);
end
fprintf(fid,'# Unit of covariance is 1e-8(rad/myr)^2 \n'); fclose(fid);

end

function [a,b,c,d,e,f]=out_cov(covmatrix)
covpol=covmatrix.*1e12.*1e8;
a=covpol(1,1);
b=covpol(1,2);
c=covpol(1,3);
d=covpol(2,2);
e=covpol(2,3);
f=covpol(3,3);
end

%% Calc optimum value from tcha.mat
function [bslip,bslipl,vec] = CalcOptimumValue(prm,tcha,G,d)
mpmean = tcha.avepol;
mcmean = tcha.aveflt;
mamean = tcha.aveasp;
mimean = tcha.aveine;

Hu = Heaviside(G(1).zulim - G(1).zc);
Hd = Heaviside(G(1).zdlim - G(1).zc);
Hlim = Hd - Hu;

idl1 = (Heaviside(G(1).zd*mamean-G(1).zc) - Heaviside(G(1).zu*mamean-G(1).zc)) .* Hlim;
idl = logical(d(1).maid *  idl1);
idc = logical(d(1).maid * ~idl1);

% Calculate back-slip on locked and creeping patches.
bslip      = (G(1).tb_mec * mpmean) .* d(1).cfinv_mec .* idl;
bslipl     = bslip;
bslip(idc) = -G(1).s(idc,idc) \ (G(1).s(idc,idl) * bslip(idl));

% Calc vectors for mean parameters
vec.rig = G(1).p * mpmean;
vec.kin = G(1).c_kin * ((G(1).tb_kin * mpmean) .* d(1).cfinv_kin .* (d(1).mcid * mcmean));
vec.mec = G(1).c_mec * bslip;
vec.ine = G(1).i * mimean;
% Zero padding
if prm.gpu ~= 99
  if isempty(vec.rig); vec.rig = zeros(size(d(1).ind),precision,'gpuArray'); end
  if isempty(vec.kin); vec.kin = zeros(size(d(1).ind),precision,'gpuArray'); end
  if isempty(vec.mec); vec.mec = zeros(size(d(1).ind),precision,'gpuArray'); end
  if isempty(vec.ine); vec.ine = zeros(size(d(1).ind),precision,'gpuArray'); end
else
  if isempty(vec.rig); vec.rig = zeros(size(d(1).ind)); end
  if isempty(vec.kin); vec.kin = zeros(size(d(1).ind)); end
  if isempty(vec.mec); vec.mec = zeros(size(d(1).ind)); end
  if isempty(vec.ine); vec.ine = zeros(size(d(1).ind)); end
end
% Total velocities
vec.sum = vec.rig + vec.kin + vec.mec + vec.ine;


end

%% Save slip deficit to txt files
function ExportSlipDeficit(folder,cal,blk,d)
mm3 = 1;
mm1 = 1;
dir_out = [folder,'/backslip'];
exid = exist(dir_out);
if exid~=7; mkdir(dir_out); end

for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    nf = size(blk(1).bound(nb1,nb2).blon,1);
    if nf~=0
      if blk(1).bound(nb1,nb2).flag2 == 1
        fid     = fopen([dir_out,'/bs_',num2str(nb1),'_',num2str(nb2),'.txt'],'w');
        clon    = mean(blk(1).bound(nb1,nb2).blon,2);
        clat    = mean(blk(1).bound(nb1,nb2).blat,2);
        cdep    = mean(blk(1).bound(nb1,nb2).bdep,2);
        sdr0    = sqrt(  cal.aslip(mm3     :mm3+  nf-1).^2 ...
                       + cal.aslip(mm3+  nf:mm3+2*nf-1).^2 ...
                       + cal.aslip(mm3+2*nf:mm3+3*nf-1).^2 );     % slip deficit before shear stress released
        sdr1    = sqrt(  cal.bslip(mm3     :mm3+  nf-1).^2 ...
                       + cal.bslip(mm3+  nf:mm3+2*nf-1).^2 ...
                       + cal.bslip(mm3+2*nf:mm3+3*nf-1).^2 );     % slip deficit after shear stress released
        flt_id  = mm1:mm1+nf-1;
        
        outdata = [flt_id' ...
                   blk(1).bound(nb1,nb2).blon ...
                   blk(1).bound(nb1,nb2).blat ...
                   blk(1).bound(nb1,nb2).bdep ...
                   clon clat cdep ...
                   sdr1 sdr0];
        
        fprintf(fid,'# %6s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %10s %10s \n',...
                        'tri_no','lon1','lon2','lon3','lat1','lat2','lat3','dep1','dep2','dep3',...
                        'c_lon','c_lat','c_dep','sdr[mm/yr]','sdr_int');
        fprintf(fid,'%8d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %10.4f %10.4f \n',outdata');
        fclose(fid);
        mm1 = mm1 +   nf;
        mm3 = mm3 + 3*nf;
      end
    end
  end
end

end

%% Export asperities
function ExportAsperities(folder,blk)
dir_out = [folder,'/backslip'];
exid = exist(dir_out);
if exid~=7; mkdir(dir_out); end
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    nf = size(blk(1).bound(nb1,nb2).blon,1);
    if nf ~= 0
      if blk(1).bound(nb1,nb2).flag2 == 1
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

end

%% Export vectors of calculation, observation, residual vector
function ExportVectors(folder,obs,cal)
% 
% Calculate vector based on sampled parameter and green function
calvec.rig = cal.rig;
calvec.ela = cal.ela;
calvec.ine = cal.ine;
calvec.sum = cal.smp;
% Calculate residual vectors
obsv = [     obs(1).evec ;      obs(1).nvec ;      obs(1).hvec ];
calv = [cal.smp(1:3:end)'; cal.smp(2:3:end)'; cal.smp(3:3:end)'];
resvec.sum=obsv-calv;
% 
% Save to txt files
odir = [folder,'/vector'];
exid = exist(odir);
if exid~=7; mkdir(odir); end
obsvec = [obs(1).alon;obs(1).alat;obs(1).evec;obs(1).nvec;obs(1).hvec];
calvec.rig = [obs(1).alon;obs(1).alat;calvec.rig(1:3:end)';calvec.rig(2:3:end)';calvec.rig(3:3:end)'];
calvec.ela = [obs(1).alon;obs(1).alat;calvec.ela(1:3:end)';calvec.ela(2:3:end)';calvec.ela(3:3:end)'];
calvec.sum = [obs(1).alon;obs(1).alat;calvec.sum(1:3:end)';calvec.sum(2:3:end)';calvec.sum(3:3:end)'];
resvec.sum = [obs(1).alon;obs(1).alat;resvec.sum];
fid = fopen([odir,'/obs_vector.txt'],'w');
fprintf(fid,'%f %f %f %f %f\n',obsvec);
fclose(fid);
fid = fopen([odir,'/cal_vector.txt'],'w');
fprintf(fid,'%f %f %f %f %f\n',calvec.sum);
fclose(fid);
fid = fopen([odir,'/cal_vector_rig_site.txt'],'w');
fprintf(fid,'%f %f %f %f %f\n',calvec.rig);
fclose(fid);
fid = fopen([odir,'/cal_vector_ela_site.txt'],'w');
fprintf(fid,'%f %f %f %f %f\n',calvec.ela);
fclose(fid);
fid = fopen([odir,'/res_vector.txt'],'w');
fprintf(fid,'%f %f %f %f %f\n',resvec.sum);
fclose(fid);
% 
end

%% Export stress and strain
function ExportStressStrain(result_dir,blk,d,G,cal)
% Stein and Wysession (2003)
lambda = 30;  % [GPa]
mu     = 30;  % [GPa]

% substitude index of locking patches
idl = logical(d(1).maid *  d(1).idl);
idc = logical(d(1).maid * ~d(1).idl);

% Calculate shear strain on meshes.
strain0 = (G(1).s(idc,idl) * cal.aslip(idl));         % pre-released shear strain
strain1 =               zeros(size(strain0));         % post-released shear strain

mm1 = 1;
mm3 = 1;
odir = [result_dir,'/stress'];
exid = exist(odir);
if exid~=7; mkdir(odir); end
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    nf = size(blk(1).bound(nb1,nb2).blon,1);
    if nf ~= 0
      if blk(1).bound(nb1,nb2).flag2 == 1
        flt_id  = mm1:mm1+nf-1;
        clon    = mean(blk(1).bound(nb1,nb2).blon,2);
        clat    = mean(blk(1).bound(nb1,nb2).blat,2);
        cdep    = mean(blk(1).bound(nb1,nb2).bdep,2);
        e.xx = 1e-6 .* zeros(nf,1);
        e.xy = 1e-6 .* zeros(nf,1);
        e.yy = 1e-6 .* zeros(nf,1);
        e.xz = 1e-6 .* strain1(mm3     :mm3+  nf-1);
        e.yz = 1e-6 .* strain1(mm3+  nf:mm3+2*nf-1);
        e.zz = 1e-6 .* strain1(mm3+2*nf:mm3+3*nf-1);
        e_int.xx = 1e-6 .* zeros(nf,1);
        e_int.xy = 1e-6 .* zeros(nf,1);
        e_int.yy = 1e-6 .* zeros(nf,1);
        e_int.xz = 1e-6 .* strain0(mm3     :mm3+  nf-1);
        e_int.yz = 1e-6 .* strain0(mm3+  nf:mm3+2*nf-1);
        e_int.zz = 1e-6 .* strain0(mm3+2*nf:mm3+3*nf-1);
        [s]     = StrainToStress(e,lambda,mu);
        [s_int] = StrainToStress(e_int,lambda,mu);
        
        fid = fopen([odir,'/stress_',num2str(nb1),'_',num2str(nb2),'.txt'],'w');
        outdata = [flt_id' ...
            blk(1).bound(nb1,nb2).blon ...
            blk(1).bound(nb1,nb2).blat ...
            blk(1).bound(nb1,nb2).bdep ...
            clon clat cdep ...
            e.xz e.yz e.zz ...
            e_int.xz e_int.yz e_int.zz ...
                s.xz     s.yz     s.zz ...
            s_int.xz s_int.yz s_int.zz];
        fprintf(fid,'# %6s %7s %7s %7s %7s %7s %7s %7s %7s %7s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s \n', ...
                    'tri_no','lon1','lon2','lon3','lat1','lat2','lat3','dep1','dep2','dep3',...
                    'c_lon','c_lat','c_dep', ...
                    'e_st','e_dp','e_ts','e_ini_st','e_ini_dp','e_ini_ts', ...
                    's_st','s_dp','s_ts','s_ini_st','s_ini_dp','s_ini_ts');
        fprintf(fid,'%8d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %10.4f %10.4f %10.4f %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e \n',outdata');
        fclose(fid);

        mm1 = mm1 +   nf;
        mm3 = mm3 + 3*nf;
      end
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

%% Heaviside step function
function y = Heaviside(x)
% Calculate Heveaside step function
%   y = Heaviside(x)
%   y and x are possible to be either scaler or vector.
y = 1 .* (x >= 0);
end

%%
function [lat,lon,ang]=xyzp2lla(X,Y,Z)
% XYZP2LLA  Converts Shpear coordinates from cartesian. Vectorized.
% GRS80
% CODE BY T.ITO 2017/03/11     ver0.1
% lat: deg, lon: deg, ang: deg/m.y.
lat=atan2(Z,sqrt(X.*X+Y.*Y)).*180/pi;
lon=atan2(Y,X).*180/pi;
ang=sqrt(X.*X+Y.*Y+Z.*Z).*(1e6.*(180./pi));
end
