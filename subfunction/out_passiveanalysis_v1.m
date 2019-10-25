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

ExportBackSlip(result_dir,cal,blk);
ExportAsperities(folder,blk);
ExportVectors(result_dir,obs,cal);
% ExportStressStrain(result_dir,blk,d,G,cal);

fprintf('Files exported. \n');
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

%% Export vectors of calculation, observation, residual vector
function ExportVectors(folder,obs,cal)
% 
% Calculate vector based on sampled parameter and green function
cal.res = reshape([obs(1).evec;obs(1).nvec;obs(1).hvec],3*obs(1).nobs,1) - cal.smp;
% 
% Save to txt files
odir = [folder,'/vector'];
exid = exist(odir);
if exid~=7; mkdir(odir); end
% Save vobs
fobs = fopen(fullfile(odir,'obs_vector.txt'),'wt');
fprintf(fobs,'# site lon  lat  hei  ve   vn   vu   se   sn   su\n');
for nob = 1:obs(1).nobs
  fprintf(fobs,'%s %f %f %f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n',...
      obs(1).name{nob},obs(1).alon(nob),obs(1).alat(nob),obs(1).ahig(nob),...
      obs(1).evec(nob),obs(1).nvec(nob),obs(1).hvec(nob),obs(1).eerr(nob),obs(1).nerr(nob),obs(1).herr(nob));
end
fclose(fobs);
% 
% Save vcal
file = fullfile(odir,'cal_vector.txt');
savevector(file,obs,cal.smp);
% Save vres
file = fullfile(odir,'res_vector.txt');
savevector(file,obs,cal.res);
% Save vrig
file = fullfile(odir,'rig_vector.txt');
savevector(file,obs,cal.rig);
% Save vmec
file = fullfile(odir,'mec_vector.txt');
savevector(file,obs,cal.mec);
% Save vkin
file = fullfile(odir,'kin_vector.txt');
savevector(file,obs,cal.kin);
end

function savevector(file,obs,v)
ve = v(1:3:end);
vn = v(2:3:end);
vu = v(3:3:end);
fid = fopen(file,'wt');
fprintf(fid,'# site lon  lat  hei  ve   vn   vu\n');
for nob = 1:obs(1).nobs
  fprintf(fid,'%s %f %f %f %6.2f %6.2f %6.2f\n',...
      obs(1).name{nob},obs(1).alon(nob),obs(1).alat(nob),obs(1).ahig(nob),...
      ve(nob),vn(nob),vu(nob));
end
fclose(fid);
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

%% Save slip deficit to txt files
function ExportBackSlip(folder,cal,blk)
savedir = fullfile(folder,'backslip');
if exist(savedir) ~=7; mkdir(savedir); end
mm1m = 1;
mm3m = 1;
mm1k = 1;
mm3k = 1;
mm1  = 1;
mm3  = 1;
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    nf = size(blk(1).bound(nb1,nb2).blon,1);
    if nf ~= 0
      fltnum = mm1:mm1+nf-1;
      clon = mean(blk(1).bound(nb1,nb2).blon,2);
      clat = mean(blk(1).bound(nb1,nb2).blat,2);
      cdep = mean(blk(1).bound(nb1,nb2).bdep,2);
      if blk(1).bound(nb1,nb2).flag2 == 1
        file = fullfile(savedir,['trimec_',num2str(nb1),'_',num2str(nb2),'.txt']);
        fmec = fopen(file,'wt');
        bslip_st  = cal.bslip( mm3m     :mm3m+  nf-1);
        bslip_dp  = cal.bslip( mm3m+  nf:mm3m+2*nf-1);
        bslip_ts  = cal.bslip( mm3m+2*nf:mm3m+3*nf-1);
        bslipl_st = cal.aslip(mm3m     :mm3m+  nf-1);
        bslipl_dp = cal.aslip(mm3m+  nf:mm3m+2*nf-1);
        bslipl_ts = cal.aslip(mm3m+2*nf:mm3m+3*nf-1);
        outdata = [fltnum',...
            blk(1).bound(nb1,nb2).blon,...
            blk(1).bound(nb1,nb2).blat,...
            blk(1).bound(nb1,nb2).bdep,...
            clon,clat,cdep,...
            bslip_st, bslip_dp, bslip_ts,...
            bslipl_st,bslipl_dp,bslipl_ts];
        fprintf(fmec,'# tri lon1 lon2 lon3 lat1 lat2 lat3 dep1 dep2 dep3 clon clat cdep bslip_st bslip_dp bslip_ts bslipl_st bslipl_dp bslipl_ts\n');
        fprintf(fmec,'%6i %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n',outdata');
        fclose(fmec);
        mm1m = mm1m +   nf;
        mm3m = mm3m + 3*nf;
      else
        %         file = fullfile(savedir,['trikin_',num2str(nb1),'_',num2str(nb2),'.txt']);
        %         fkin = fopen(file,'wt');
        %         sdr_st  = bslip.kin(mm3k     :mm3k+  nf-1);
        %         sdr_dp  = bslip.kin(mm3k+  nf:mm3k+2*nf-1);
        %         sdr_ts  = bslip.kin(mm3k+2*nf:mm3k+3*nf-1);
        %         cpmean  = tcha.aveflt(mm1k:mm1k+nf-1);
        %         outdata = [fltnum',...
        %             blk(1).bound(nb1,nb2).blon,...
        %             blk(1).bound(nb1,nb2).blat,...
        %             blk(1).bound(nb1,nb2).bdep,...
        %             clon,clat,cdep,...
        %             cpmean, sdr_st, sdr_dp, sdr_ts];
        %         fprintf(fkin,'# tri lon1 lon2 lon3 lat1 lat2 lat3 dep1 dep2 dep3 clon clat cdep coupling sdr_st sdr_dp sdr_ts\n');
        %         fprintf(fkin,'%6i %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %10.4f %10.4f %10.4f %10.4f\n',outdata');
        %         fclose(fkin);
        mm1k = mm1k +   nf;
        mm3k = mm3k + 3*nf;
      end
      mm1 = mm1 +   nf;
      mm3 = mm3 + 3*nf;
    end
  end
end

end