%% Main part (single MCMC results)
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
% Common export
[blk] = ReadAsperityRegions(savedir,blk,prm);
[blk] = AsperityPoint(blk,obs);
SaveBlockLine(savedir,blk)
SaveAsperitySegmentArea(savedir,blk,obs,tcha)
SaveAsperityPoint(savedir,blk);
% Export for each method
for rep = 1:prm.nrep
  if strcmpi(prm.method,'MCMC.RE')
    savedir = fullfile(pwd,folder,['T_',num2str(rep,'%02i')]);
  end
  [bslip,slip,vec] = CalcOptimumValue(prm,obs,tcha,G,d,rep);
  SavePoles(savedir,blk,tcha,rep);
  SaveBackslip(savedir,blk,tcha,bslip,slip,rep);
  SaveVectors(savedir,obs,vec);
  SaveInternalStrain(savedir,tcha,blk,prm,obs,G,rep);
  SaveRelativeMotion(savedir,blk,tcha,rep);
end

fprintf('Files exported. \n');
end

%% Save relative motion
function SaveRelativeMotion(folder,blk,tcha,rep)
folder=[folder,'/vector'];
if exist(folder)~=7; mkdir(folder); end
fid=fopen([folder,'/relative_motion.txt'],'w');
fprintf(fid,'# %5s %7s %7s %7s %10s %10s %10s \n',...
    'Lon1','Lon2','Lat1','Lat2','abs_Vel','str_Vel','dip_Vel');
for nb1=1:blk(1).nblock
  blk(nb1).pol=[tcha.avepol(3.*nb1-2,:,rep);tcha.avepol(3.*nb1-1,:,rep);tcha.avepol(3.*nb1,:,rep)];
  for nb2=nb1+1:blk(1).nblock
    blk(nb2).pol(:)=[tcha.avepol(3.*nb2-2,:,rep);tcha.avepol(3.*nb2-1,:,rep);tcha.avepol(3.*nb2,:,rep)];
    if ~isempty(blk(1).bound(nb1,nb2).lat)
      fprintf(fid,'> %s - %s \n',num2str(nb1),num2str(nb2));
      vel = CalcRelativeMotion(blk,nb1,nb2);
      fprintf(fid,'%7.3f %7.3f %7.3f %7.3f %10.4f %10.4f %10.4f \n',...
          [blk(1).bound(nb1,nb2).lon(1:end-1), blk(1).bound(nb1,nb2).lon(2:end),...
           blk(1).bound(nb1,nb2).lat(1:end-1), blk(1).bound(nb1,nb2).lat(2:end),...
           vel.abs, vel.str, vel.dip]');
    end
  end
end
% 
end

function [vel] = CalcRelativeMotion(blk,nb1,nb2)
blon0 = mean(blk(1).bound(nb1,nb2).lon);
blat0 = mean(blk(1).bound(nb1,nb2).lat);
[bxy(:,1),bxy(:,2)] = PLTXY(blk(nb1).lat,blk(nb1).lon,blat0,blon0);
[boxy(:,1),boxy(:,2)] = PLTXY(blk(1).bound(nb1,nb2).lat,blk(1).bound(nb1,nb2).lon,blat0,blon0);
bocxy = (boxy(1:end-1,:) + boxy(2:end,:)) ./ 2;
boclat = (blk(1).bound(nb1,nb2).lat(1:end-1)+blk(1).bound(nb1,nb2).lat(2:end))./2;
boclon = (blk(1).bound(nb1,nb2).lon(1:end-1)+blk(1).bound(nb1,nb2).lon(2:end))./2;
bocxyz = conv2ell(boclat,boclon);
uv = [0 0 1];
boline = boxy(2:end,:)-boxy(1:end-1,:);
bonormal = [uv(2).*          0 - uv(3).*boline(:,2),...
            uv(3).*boline(:,1) - uv(1).*0          ,...
            uv(1).*boline(:,2) - uv(2).*boline(:,1)];
v12 = pole2velo((blk(nb2).pol(:)-blk(nb1).pol(:))',bocxyz);
v12x = v12(1:2:end);
v12y = v12(2:2:end);
v12str = (v12x.*boline(:,1)+v12y.*boline(:,2)).*boline...
           ./(boline(:,1).^2+boline(:,2).^2);
v12dip = (v12x.*bonormal(:,1)+v12y.*bonormal(:,2)).*bonormal(:,1:2)...
           ./(bonormal(:,1).^2+bonormal(:,2).^2);
v12st = zeros(size(v12x));
v12dp = zeros(size(v12x));
v12abs = sqrt(v12x.^2+v12y.^2);
for nbo = 1:size(blk(1).bound(nb1,nb2).lat,1)-1
  st = [v12str(nbo,:) 0]; v12st(nbo) = sqrt(st(1).^2+st(2).^2);
  dp = [v12dip(nbo,:) 0]; v12dp(nbo) = sqrt(dp(1).^2+dp(2).^2);
  nv = cross(uv,st)./norm(cross(uv,st),2);
  dp = dp./norm(dp,2);
  lstr = bocxy(nbo,:)+nv(1:2);
  ldip = bocxy(nbo,:)+dp(1:2);
  inidst = inpolygon(lstr(1),lstr(2),bxy(:,1),bxy(:,2));
  iniddp = inpolygon(ldip(1),ldip(2),bxy(:,1),bxy(:,2));
  if inidst == 1; v12st(nbo) = -1*v12st(nbo); end % Left lateral
  if iniddp ~= 1; v12dp(nbo) = -1*v12dp(nbo); end % Open
end
vel.abs = v12abs;
vel.str = v12st;
vel.dip = v12dp;
end

function Vneu = pole2velo(Pxyz,Oxyz)
% pole2velo Convert velocity from Euler pole. Vectorized.
[Nobs,~]=size(Oxyz);
[Npol,~]=size(Pxyz);
Vxyz=zeros(Npol,3,'single'); 
Vneu=zeros(2.*Nobs,Npol,'single');
%
for N=1:Nobs
  Vxyz(:,1) = -Oxyz(N,2).*Pxyz(:,3) + Pxyz(:,2).*Oxyz(N,3);
  Vxyz(:,2) = -Oxyz(N,3).*Pxyz(:,1) + Pxyz(:,3).*Oxyz(N,1);
  Vxyz(:,3) = -Oxyz(N,1).*Pxyz(:,2) + Pxyz(:,1).*Oxyz(N,2);
  Vneu(2.*N-1,:) =            -Oxyz(N,5).*Vxyz(:,1) ...
                              +Oxyz(N,7).*Vxyz(:,2); %E
  Vneu(2.*N,:)   = -Oxyz(N,4).*Oxyz(N,7).*Vxyz(:,1) ...
                   -Oxyz(N,4).*Oxyz(N,5).*Vxyz(:,2) ...
                              +Oxyz(N,6).*Vxyz(:,3); %N
end
end

%% Save internal strains
function SaveInternalStrain(folder,tcha,blk,prm,obs,G,rep)
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
vinternal = G(1).i * tcha.aveine(:,:,rep);
outdata = [obs(1).alat; obs(1).alon; vinternal(1:3:end)'; vinternal(2:3:end)'];
fprintf(fidv,'%7.3f %7.3f %10.4f %10.4f \n',outdata); fclose(fidv);
clear vinternal outdata
for nb = 1:blk(1).nblock
  Gb.i = zeros(size(G.i));
  Gb.i(:,3*nb-2:3*nb) = G.i(:,3*nb-2:3*nb);
  vinternal = Gb.i * tcha.aveine(:,:,rep);  % Displacement due to internal deformation
  outdata = [obs(1).alat; obs(1).alon; vinternal(1:3:end)'; vinternal(2:3:end)'];
  fidv = fopen([savefolder,'/Internal_Deformation_vector_blk',num2str(nb),'.txt'],'w');
  fprintf(fidv,'%7.3f %7.3f %10.4f %10.4f \n',outdata);
  fclose(fidv);
  exx = tcha.aveine(3*nb-2,:,rep);
  exy = tcha.aveine(3*nb-1,:,rep);
  eyy = tcha.aveine(3*nb  ,:,rep);
  sigexx = tcha.stdine(3*nb-2,:,rep);
  sigexy = tcha.stdine(3*nb-1,:,rep);
  sigeyy = tcha.stdine(3*nb  ,:,rep);
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
% Save vrig
file = fullfile(savedir,'rig_vector.txt');
savevector(file,obs,vec.rig);
% Save vmec
file = fullfile(savedir,'mec_vector.txt');
savevector(file,obs,vec.mec);
% Save vkin
file = fullfile(savedir,'kin_vector.txt');
savevector(file,obs,vec.kin);
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
function SaveBackslip(folder,blk,tcha,bslip,slip,rep)
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
        slip_st   =  slip.mec( mm3m     :mm3m+  nf-1);
        slip_dp   =  slip.mec( mm3m+  nf:mm3m+2*nf-1);
        slip_ts   =  slip.mec( mm3m+2*nf:mm3m+3*nf-1);
        slipl_st  =  slip.mecl(mm3m     :mm3m+  nf-1);
        slipl_dp  =  slip.mecl(mm3m+  nf:mm3m+2*nf-1);
        slipl_ts  =  slip.mecl(mm3m+2*nf:mm3m+3*nf-1);
        pc = tcha.aveaid(mm1m:mm1m+nf-1,:,rep);
        outdata = [fltnum',...
            blk(1).bound(nb1,nb2).blon,...
            blk(1).bound(nb1,nb2).blat,...
            blk(1).bound(nb1,nb2).bdep,...
            clon,clat,cdep,...
            bslip_st, bslip_dp, bslip_ts,...
            bslipl_st,bslipl_dp,bslipl_ts,...
            slip_st,slip_dp,slip_ts,...
            slipl_st,slipl_dp,slipl_ts,...
            pc];
        fprintf(fmec,'#   1    2    3    4    5    6    7    8    9   10   11   12   13       14       15       16        17        18        19      20      21      22       23       24       25         26\n');
        fprintf(fmec,'# tri lon1 lon2 lon3 lat1 lat2 lat3 dep1 dep2 dep3 clon clat cdep bslip_st bslip_dp bslip_ts bslipl_st bslipl_dp bslipl_ts slip_st slip_dp slip_ts slipl_st slipl_dp slipl_ts p_coupling\n');
        fprintf(fmec,'%6i %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n',outdata');
        fclose(fmec);
        mm1m = mm1m +   nf;
        mm3m = mm3m + 3*nf;
      else
        file = fullfile(savedir,['trikin_',num2str(nb1),'_',num2str(nb2),'.txt']);
        fkin = fopen(file,'wt');
        sdr_st  = bslip.kin(mm3k     :mm3k+  nf-1);
        sdr_dp  = bslip.kin(mm3k+  nf:mm3k+2*nf-1);
        sdr_ts  = bslip.kin(mm3k+2*nf:mm3k+3*nf-1);
        sr_st   =  slip.kin(mm3k     :mm3k+  nf-1);
        sr_dp   =  slip.kin(mm3k+  nf:mm3k+2*nf-1);
        sr_ts   =  slip.kin(mm3k+2*nf:mm3k+3*nf-1);
        cpmean  = tcha.aveflt(mm1k:mm1k+nf-1,:,rep);
        outdata = [fltnum',...
            blk(1).bound(nb1,nb2).blon,...
            blk(1).bound(nb1,nb2).blat,...
            blk(1).bound(nb1,nb2).bdep,...
            clon,clat,cdep,...
            cpmean, sdr_st, sdr_dp, sdr_ts,...
            sr_st, sr_dp, sr_ts];
        fprintf(fkin,'#   1    2    3    4    5    6    7    8    9   10   11   12   13       14     15     16     17    18    19    20\n');
        fprintf(fkin,'# tri lon1 lon2 lon3 lat1 lat2 lat3 dep1 dep2 dep3 clon clat cdep coupling sdr_st sdr_dp sdr_ts sr_st sr_dp sr_ts\n');
        fprintf(fkin,'%6i %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n',outdata');
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
function SavePoles(folder,blk,tcha,rep)
savedir = fullfile(folder,'block');
if exist(savedir) ~= 7; mkdir(savedir); end
fid = fopen(fullfile(savedir,'/est_euler_pole.txt'),'wt');
fprintf(fid,'# %6s %8s %7s %8s %9s %9s %9s %9s %9s %9s \n',...
    'Blk_no','Lat(deg)','Lon(deg)','Ang(deg/my)',...
    'sigxx','sigxy','sigxz','sigyy','sigyz','sigzz');
fprintf('# %6s %8s %7s %8s %9s %9s %9s %9s %9s %9s \n',...
    'Blk_no','Lat(deg)','Lon(deg)','Ang(deg/my)',...
    'sigxx','sigxy','sigxz','sigyy','sigyz','sigzz');
for bk = 1:blk(1).nblock
  [latp,lonp,ang] = xyzp2lla(tcha.avepol(3.*bk-2,:,rep),tcha.avepol(3.*bk-1,:,rep),tcha.avepol(3.*bk,:,rep));
  [a,b,c,d,e,f] = out_cov(tcha.covpol(3*bk-2:3*bk,3*bk-2:3*bk,rep));
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

%% calc asperity segment area
function SaveAsperitySegmentArea(folder,blk,obs,tcha)
savedir = fullfile(folder,'backslip');
if exist(savedir) ~=7; mkdir(savedir); end
alat0 = mean(obs(1).alat,2);
alon0 = mean(obs(1).alon,2);
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
      [trix,triy] = PLTXY(blk(1).bound(nb1,nb2).blat,blk(1).bound(nb1,nb2).blon,alat0,alon0);
      triz = blk(1).bound(nb1,nb2).bdep;
      [tricx,tricy] = PLTXY(mean(blk(1).bound(nb1,nb2).blat,2),mean(blk(1).bound(nb1,nb2).blon,2),alat0,alon0);
      tricz = mean(blk(1).bound(nb1,nb2).bdep,2);
      area = zeros(size(blk(1).bound(nb1,nb2).blat,1),1);
      for ntri = 1:size(blk(1).bound(nb1,nb2).blat,1)
        area(ntri) = triangle_area([trix(ntri,:)', triy(ntri,:)', triz(ntri,:)']);
      end
      blk(1).bound(nb1,nb2).triarea = area;
      if blk(1).bound(nb1,nb2).flag2 == 1
        if blk(1).bound(nb1,nb2).segmentid == 1
          for nseg = 1:size(blk(1).bound(nb1,nb2).segment,2)
            [edgex,edgey] = PLTXY(blk(1).bound(nb1,nb2).segment(nseg).lat,blk(1).bound(nb1,nb2).segment(nseg).lon,alat0,alon0);
            segid = inpolygon(tricx,tricy,edgex,edgey);
            asp(1).bound(nb1,nb2).segment(nseg).smparea = segid'.*blk(1).bound(nb1,nb2).triarea' * tcha.smpaid(mm1m:mm1m+nf-1,:);
          end
        else
          asp(1).bound(nb1,nb2).smparea = blk(1).bound(nb1,nb2).triarea' * tcha.smpaid(mm1m:mm1m+nf-1,:);
        end
        mm1m = mm1m +   nf;
        mm3m = mm3m + 3*nf;  
      else
        mm1k = mm1k +   nf;
        mm3k = mm3k + 3*nf;          
      end
      mm1 = mm1 +   nf;
      mm3 = mm3 + 3*nf;
    end
  end
end
save(fullfile(savedir,'asp'),'asp');
end

%% Save block lines
function SaveBlockLine(folder,blk)
savedir = fullfile(folder,'block');
if exist(savedir) ~= 7; mkdir(savedir); end
fid = fopen(fullfile(savedir,'BlockLine.txt'),'wt');
for nb = 1:size(blk,2)
  blkname = strsplit(blk(nb).name,'.txt');
  fprintf(fid,'> %s\n',char(blkname{1}));
  fprintf(fid,'%f %f\n',[blk(nb).lon,blk(nb).lat]');
end
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
function [Bslip,Slip,vec] = CalcOptimumValue(prm,obs,tcha,G,d,rep)
mpmean = tcha.avepol(:,:,rep);
mcmean = tcha.aveflt(:,:,rep);
mamean = tcha.aveasp(:,:,rep);
mimean = tcha.aveine(:,:,rep);

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
Slip.mec   = -((G(1).tb_mec * mpmean) .* d(1).cfinv_mec - Bslip.mec );
Slip.mecl  = -((G(1).tb_mec * mpmean) .* d(1).cfinv_mec - Bslip.mecl); Slip.mecl(idc) = NaN;
Slip.kin   = -((G(1).tb_kin * mpmean) .* d(1).cfinv_kin - Bslip.kin );

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

%% Read asperity regions
function [blk] = ReadAsperityRegions(savedir,blk,prm)
%    Test version coded by H. Kimura 2019/1/29
% Revised version coded by H. Kimura 2019/2/5

for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    if blk(1).bound(nb1,nb2).flag2 == 1
      blk(1).bound(nb1,nb2).segmentid = 0;
      patchfile = fullfile(savedir,'backslip',['segments_',num2str(nb1),'_',num2str(nb2),'.txt']);
      fid       = fopen(patchfile,'r');
      if fid >= 0
        blk(1).bound(nb1,nb2).segmentid = 1;
        np    = 0;
        n     = 0;
        while 1
          tline = fgetl(fid);
          if ~ischar(tline) ; break; end
          if tline(1) ~= '>'
            n   = n+1;
            tmp = strsplit(strtrim(tline));
            blk(1).bound(nb1,nb2).segment(np+1).lon(n) = str2double(cellstr(tmp(1)));
            blk(1).bound(nb1,nb2).segment(np+1).lat(n) = str2double(cellstr(tmp(2)));
          else
            np = np+1;
            n  = 0;
            continue;
          end
        end
      else
        %         error(['Not found', patchfile]);
      end
    end
  end
end

fprintf('=== Read Locked Patches=== \n');
end

%%
function [area]=triangle_area(P,method)
% This function gives the area of a triangle
%
% [area]=triangle_area(Points, Method)
%
% Points: The Points should be a numeric array, of size 3xn, 
%         thus the points can be 2D, 3D... nD
% Method: Can be 'h' area calculation with Heron's formula 
%         or can be 'q' Orthogonal-triangular decomposition (default)
%
% Example: 
% P1=[0 0]; P2=[1 0.5]; P3=[0.5 1];
% area = triangle_area([P1;P2;P3])
%
% Version 1.1 updated on 2007-09-21 
% Added 'Orthogonal-triangular decomposition' after a usefull review of John D'Errico 

% Default output format
if(exist('method','var')==0), method='q'; end

% Check input
if((method~='h')&&(method~='q')), error('Unknown area calculation method'); end
[k,m]=size(P); if(k~=3), error('Points are not a 3xn array'); end

if(method=='h')
    % Length of edges
    L=[sqrt(sum((P(1,:)-P(2,:)).^2)) sqrt(sum((P(2,:)-P(3,:)).^2)) sqrt(sum((P(3,:)-P(1,:)).^2))];
    
    % Area calculation with Heron's formula
    s = ((L(1)+L(2)+L(3))/2); 
    area = sqrt(s*(s-L(1))*(s-L(2))*(s-L(3)));
else
    % Area calculation with Orthogonal-triangular decomposition
    [q,r] = qr((P(2:3,:) - repmat(P(1,:),2,1))');
    area=abs(prod(diag(r)))/2;
end
    
end

%% Heaviside step function
function y = Heaviside(x)
% Calculate Heveaside step function
%   y = Heaviside(x)
%   y and x are possible to be either scaler or vector.
y = 1 .* (x >= 0);
end

%% Coordinates conversion
function [x,y,z]=ell2xyz(lat,lon,h)
% ELL2XYZ  Converts ellipsoidal coordinates to cartesian. Vectorized.
% GRS80
% CODE BY T.ITO 2006/12/13     ver0.1
% BUG FIX  BY T.ITO 2015/11/13 ver0.2
% 
a=6378137.0; % m
f=1./298.257222101;
e2=1-(1-f)^2;
%
rad=pi/180;
lat=lat.*rad;
lon=lon.*rad;
clat=cos(lat);
clon=cos(lon);
slat=sin(lat);
slon=sin(lon);
%
v=a./sqrt(1-e2.*slat.*slat);
x=(v+h).*clat.*clon;
y=(v+h).*clat.*slon;
z=(v.*(1-e2)+h).*slat;
end

function [OOxyz]=conv2ell(Olat,Olon)
Olat=Olat(:);
Olon=Olon(:);
deg2rad=pi/180;
[Oxyz(:,1),Oxyz(:,2),Oxyz(:,3)]=ell2xyz(Olat,Olon,0);
Oxyz = Oxyz*1e3;
OOxyz=[Oxyz sin(Olat*deg2rad) sin(Olon*deg2rad) cos(Olat*deg2rad) cos(Olon*deg2rad)];
end

function [lat,lon,ang]=xyzp2lla(X,Y,Z)
% XYZP2LLA  Converts Shpear coordinates from cartesian. Vectorized.
% GRS80
% CODE BY T.ITO 2017/03/11     ver0.1
% lat: deg, lon: deg, ang: deg/m.y.
lat=atan2(Z,sqrt(X.*X+Y.*Y)).*180/pi;
lon=atan2(Y,X).*180/pi;
ang=sqrt(X.*X+Y.*Y+Z.*Z).*(1e6.*(180./pi));
end

function [X,Y]=PLTXY(ALAT,ALON,ALAT0,ALON0)
%-------------------
%  PLTXY TRANSFORMS (ALAT,ALONG) TO (X,Y)
%  WHEN ICORD.NE.0  PLTXY MAKES NO CHANGE IN 
%  TRANSFORMATION BETWEEN (X,Y) AND (ALAT,ALONG).
%-------------------
A=6.378160e3;
E2=6.6944541e-3;
E12=6.7395719e-3;
D=5.72958e1;
RD=1.0/D;
RLAT = RD.*ALAT;
SLAT = sin(RLAT);
CLAT = cos(RLAT);
V2   = 1.0 + E12.*CLAT.^2;
AL   = ALON-ALON0;
PH1  = ALAT + (V2.*AL.^2.*SLAT.*CLAT)./(2.0*D);
RPH1 = PH1.*RD;
RPH2 = (PH1 + ALAT0).*0.5.*RD;
R    = A.*(1.0-E2)./sqrt((1.0-E2.*sin(RPH2).^2).^3);
AN   = A./sqrt(1.0-E2.*sin(RPH1).^2);
C1   = D./R;
C2   = D./AN;
Y    = (PH1-ALAT0)./C1;
X    = (AL.*CLAT)./C2+(AL.^3.*CLAT.*cos(2.0.*RLAT))./(6.0.*C2.*D.^2);
end
