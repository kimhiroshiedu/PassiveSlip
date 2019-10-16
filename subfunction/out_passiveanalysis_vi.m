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

[bslip,vec] = CalcOptimumValue(prm,obs,tcha,G,d);
[blk] = AsperityPoint(blk,obs);
SaveAsperityPoint(savedir,blk)
SavePoles(savedir,blk,tcha);
SaveBackslip(savedir,blk,tcha,bslip);
SaveVectors(savedir,obs,vec);
SaveInternalStrain(savedir,tcha,blk,prm,obs,G);

fprintf('Files exported. \n');
end

%% Save internal strains
function SaveInternalStrain(folder,tcha,blk,prm,obs,G)
if isfield(tcha,'aveine')==0; return; end  % old type

NN=1;
savefolder=[folder,'/innerdeform'];
exid=exist(savefolder);
if exid~=7; mkdir(savefolder); end

fide = fopen([savefolder,'/Internal_Deformation_strain.txt'],'w');
fidv = fopen([savefolder,'/Internal_Deformation_vector.txt'],'w');
fprintf(fide,'# %5s %7s %7s %7s %7s %7s %7s %7s %10s %10s %10s %10s %10s %10s %10s %10s %s\n # Unit of strain is [nanostrain/yr] \n',...
                  'Block','Lat','Lat','exx','exy','eyy','emax','emin','thetaP','shearMAX',...
                  'sig_exx','sig_exy','sig_eyy','sig_emax','sig_emin','sig_thetaP','sig_shearMAX');
fprintf(fidv,'# Latitude Longitude VE VN \n');
vinternal = G(1).i * tcha.aveine;
outdata = [obs(1).alat; obs(1).alon; vinternal(1:3:end)'; vinternal(2:3:end)'];
fprintf(fidv,'%7.3f %7.3f %10.4f %10.4f \n',outdata); fclose(fidv);
clear vinternal outdata
for nb = 1:blk(1).nblock
  Gb.i = zeros(size(G.i));
  Gb.i(:,3*nb-2:3*nb) = G.i(:,3*nb-2:3*nb);
  vinternal = Gb.i * tcha.aveine;  % Displacement due to internal deformation
  outdata = [obs(1).alat; obs(1).alon; vinternal(1:3:end)'; vinternal(2:3:end)'];
  fidv = fopen([savefolder,'/Internal_Deformation_vector_blk',num2str(nb),'.txt'],'w');
  fprintf(fidv,'%7.3f %7.3f %10.4f %10.4f \n',outdata);
  fclose(fidv);
  exx = tcha.aveine(3*nb-2);
  exy = tcha.aveine(3*nb-1);
  eyy = tcha.aveine(3*nb  );
  sigexx = tcha.stdine(3*nb-2);
  sigexy = tcha.stdine(3*nb-1);
  sigeyy = tcha.stdine(3*nb  );
  E = [exx exy;...
       exy eyy];
  [eigv,eigd] = eig(E);
  e1 = eigd(1,1); e2 = eigd(2,2);
  v1 = eigv(:,1); v2 = eigv(:,2);
  if e1 >= e2
    emax = e1; axmax = v1;
    emin = e2; axmin = v2;
  else
    emax = e2; axmax = v2;
    emin = e1; axmin = v1;
  end
  thetaP = 90-rad2deg(atan2(axmax(2),axmax(1))); % Azimth from North
  if thetaP<0; thetaP = thetaP+360; end
  shearMAX    = sqrt((1/4)*(exx-eyy)^2+exy^2);
  sigemax     = sqrt( ( 0.5 + 0.25*( (exx-eyy)^2 /4 + exy^2 )^-0.5 *( exx-eyy ) )^2 *sigexx^2 ...
                     +(          1*( (exx-eyy)^2 /4 + exy^2 )^-0.5 *  exy       )^2 *sigexy^2 ...
                     +( 0.5 - 0.25*( (exx-eyy)^2 /4 + exy^2 )^-0.5 *( exx-eyy ) )^2 *sigeyy^2 );
  sigemin     = sqrt( ( 0.5 - 0.25*( (exx-eyy)^2 /4 + exy^2 )^-0.5 *( exx-eyy ) )^2 *sigexx^2 ...
                     +(         -1*( (exx-eyy)^2 /4 + exy^2 )^-0.5 *  exy       )^2 *sigexy^2 ...
                     +( 0.5 + 0.25*( (exx-eyy)^2 /4 + exy^2 )^-0.5 *( exx-eyy ) )^2 *sigeyy^2 );
  sigshearMAX = sqrt( (       0.25*( (exx-eyy)^2 /4 + exy^2 )^-0.5 *( exx-eyy ) )^2 *sigexx^2 ...
                     +(          1*( (exx-eyy)^2 /4 + exy^2 )^-0.5 *  exy       )^2 *sigexy^2 ...
                     +(      -0.25*( (exx-eyy)^2 /4 + exy^2 )^-0.5 *( exx-eyy ) )^2 *sigeyy^2 );
  sigthetaP   = sqrt( (     -exy  / ( 4*exy^2 - (exx-eyy)^2 )                   )^2 *sigexx^2 ...
                     +( (exx-eyy) / ( 4*exy^2 - (exx-eyy)^2 )                   )^2 *sigexy^2 ...
                     +(      exy  / ( 4*exy^2 - (exx-eyy)^2 )                   )^2 *sigeyy^2 );
  fprintf(fide,'%7i %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n',...
                    nb, blk(nb).latinter, blk(nb).loninter,...
                    exx*1e9, exy*1e9, eyy*1e9,...
                    emax*1e9, emin*1e9,...
                    thetaP, shearMAX*1e9,...
                    sigexx*1e9, sigexy*1e9, sigeyy*1e9,...
                    sigemax*1e9, sigemin*1e9,...
                    rad2deg(sigthetaP), sigshearMAX*1e9);
  S(nb).exx = exx;
  S(nb).exy = exy;
  S(nb).eyy = eyy;
  S(nb).sig_exx = sigexx;
  S(nb).sig_exy = sigexy;
  S(nb).sig_eyy = sigeyy;
  S(nb).e1 = e1;
  S(nb).e2 = e2;
  S(nb).v1 = v1;
  S(nb).v2 = v2;
  S(nb).p_theta = thetaP;
  S(nb).shearMAX = shearMAX;
  S(nb).sig_emax = sigemax;
  S(nb).sig_emin = sigemin;
  S(nb).sig_thetaP = rad2deg(sigthetaP);
  S(nb).sig_shearmax = sigshearMAX;
end
fclose(fide);
strainfile=[savefolder,'/strain'];
save(strainfile,'S','-v7.3')
end

%% Save obs-, cal-, and res-vectors
function SaveVectors(folder,obs,vec)
savedir = fullfile(folder,'vector');
if exist(savedir) ~=7; mkdir(savedir); end
% Save vobs
fobs = fopen(fullfile(savedir,'obs_vector.txt'),'wt');
fprintf(fobs,'# site lon  lat  hei  ve   vn   vu   se   sn   su\n');
for nob = 1:obs(1).nobs
  fprintf(fobs,'%s %f %f %f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n',...
      obs(1).name{nob},obs(1).alon(nob),obs(1).alat(nob),obs(1).ahig(nob),...
      obs(1).evec(nob),obs(1).nvec(nob),obs(1).hvec(nob),obs(1).eerr(nob),obs(1).nerr(nob),obs(1).herr(nob));
end
fclose(fobs);
% Save vcal
file = fullfile(savedir,'cal_vector.txt');
savevector(file,obs,vec.sum);
% Save vres
file = fullfile(savedir,'res_vector.txt');
savevector(file,obs,vec.res);

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

%% Save coupling and backslip
function SaveBackslip(folder,blk,tcha,bslip)
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
        bslip_st  = bslip.mec( mm3m     :mm3m+  nf-1);
        bslip_dp  = bslip.mec( mm3m+  nf:mm3m+2*nf-1);
        bslip_ts  = bslip.mec( mm3m+2*nf:mm3m+3*nf-1);
        bslipl_st = bslip.mecl(mm3m     :mm3m+  nf-1);
        bslipl_dp = bslip.mecl(mm3m+  nf:mm3m+2*nf-1);
        bslipl_ts = bslip.mecl(mm3m+2*nf:mm3m+3*nf-1);
        pc = tcha.aveaspid(mm1m:mm1m+nf-1);
        outdata = [fltnum',...
            blk(1).bound(nb1,nb2).blon,...
            blk(1).bound(nb1,nb2).blat,...
            blk(1).bound(nb1,nb2).bdep,...
            clon,clat,cdep,...
            bslip_st, bslip_dp, bslip_ts,...
            bslipl_st,bslipl_dp,bslipl_ts,...
            pc];
        fprintf(fmec,'# tri lon1 lon2 lon3 lat1 lat2 lat3 dep1 dep2 dep3 clon clat cdep bslip_st bslip_dp bslip_ts bslipl_st bslipl_dp bslipl_ts p_coupling\n');
        fprintf(fmec,'%6i %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n',outdata');
        fclose(fmec);
        mm1m = mm1m +   nf;
        mm3m = mm3m + 3*nf;
      else
        file = fullfile(savedir,['trikin_',num2str(nb1),'_',num2str(nb2),'.txt']);
        fkin = fopen(file,'wt');
        sdr_st  = bslip.kin(mm3k     :mm3k+  nf-1);
        sdr_dp  = bslip.kin(mm3k+  nf:mm3k+2*nf-1);
        sdr_ts  = bslip.kin(mm3k+2*nf:mm3k+3*nf-1);
        cpmean  = tcha.aveflt(mm1k:mm1k+nf-1);
        outdata = [fltnum',...
            blk(1).bound(nb1,nb2).blon,...
            blk(1).bound(nb1,nb2).blat,...
            blk(1).bound(nb1,nb2).bdep,...
            clon,clat,cdep,...
            cpmean, sdr_st, sdr_dp, sdr_ts];
        fprintf(fkin,'# tri lon1 lon2 lon3 lat1 lat2 lat3 dep1 dep2 dep3 clon clat cdep coupling sdr_st sdr_dp sdr_ts\n');
        fprintf(fkin,'%6i %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %10.4f %10.4f %10.4f %10.4f\n',outdata');
        fclose(fkin);
        mm1k = mm1k +   nf;
        mm3k = mm3k + 3*nf;
      end
      mm1 = mm1 +   nf;
      mm3 = mm3 + 3*nf;
    end
  end
end
end

%% Save Euler vectors
function SavePoles(folder,blk,tcha)
fid = fopen(fullfile(folder,'/est_euler_pole.txt'),'wt');
fprintf(fid,'# %6s %8s %7s %8s %9s %9s %9s %9s %9s %9s \n',...
    'Blk_no','Lat(deg)','Lon(deg)','Ang(deg/my)',...
    'sigxx','sigxy','sigxz','sigyy','sigyz','sigzz');
fprintf('# %6s %8s %7s %8s %9s %9s %9s %9s %9s %9s \n',...
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

%% Save Asperity Points
function SaveAsperityPoint(folder,blk)
savedir = fullfile(folder,'backslip');
if exist(savedir) ~=7; mkdir(savedir); end
fid = fopen(fullfile(savedir,'/udline.txt'),'wt');
fprintf(fid,'# N      lon_d      lat_d      lon_u      lat_u \n');
lineid = [1:size(blk(1).aline_lonu,1)]';
fprintf(fid,'%3i %10.4f %10.4f %10.4f %10.4f \n',...
            [lineid,blk(1).aline_lond,blk(1).aline_latd,blk(1).aline_lonu,blk(1).aline_latu]');
fclose(fid);
end

%% Show asperity edge point indices
function [blk] = AsperityPoint(blk,obs)
blk(1).naspline  =  0;
blk(1).aline_lonu = [];
blk(1).aline_lond = [];
blk(1).aline_latu = [];
blk(1).aline_latd = [];
%
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    nf = size(blk(1).bound(nb1,nb2).blon,1);
    if nf ~= 0
      if blk(1).bound(nb1,nb2).flag2 == 1
        blk(1).aline_lonu = [blk(1).aline_lonu; blk(1).bound(nb1,nb2).asp_lonu];
        blk(1).aline_lond = [blk(1).aline_lond; blk(1).bound(nb1,nb2).asp_lond];
        blk(1).aline_latu = [blk(1).aline_latu; blk(1).bound(nb1,nb2).asp_latu];
        blk(1).aline_latd = [blk(1).aline_latd; blk(1).bound(nb1,nb2).asp_latd];
      end
    end
  end
end
fprintf('=== Asperity edge points, below === \n');
end

%% Calc optimum value from tcha.mat
function [Bslip,vec] = CalcOptimumValue(prm,obs,tcha,G,d)
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
Bslip.mec  = bslip;
Bslip.mecl = bslipl;
Bslip.kin  = ((G(1).tb_kin * mpmean) .* d(1).cfinv_kin .* (d(1).mcid * mcmean));

% Calc vectors for mean parameters
vec.rig = G(1).p * mpmean;
vec.kin = G(1).c_kin * Bslip.kin;
vec.mec = G(1).c_mec * Bslip.mec;
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
vobs = reshape([obs(1).evec; obs(1).nvec; obs(1).hvec],3*obs(1).nobs,1);
vec.res = vobs - vec.sum;

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
