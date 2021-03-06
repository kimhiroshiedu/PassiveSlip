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
% [blk] = AsperityPoint(blk,obs);
SaveBlockLine(savedir,blk)
SaveAsperitySegmentArea(savedir,blk,obs,tcha)
SaveUDlines(savedir,blk);
% Export for each method
for rep = 1:prm.nrep
  if strcmpi(prm.method,'MCMC.RE')
    savedir = fullfile(pwd,folder,['T_',num2str(rep,'%02i')]);
  end
  [s,vec] = CalcOptimumValue(prm,obs,tcha,G,d,rep);
  SavePoles(savedir,blk,tcha,rep);
  SaveBackslip(savedir,blk,tcha,s,rep);
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
savevector3(file,obs,vec.sum,vec.sum25,vec.sum50,vec.sum75);
% Save vres
file = fullfile(savedir,'res_vector.txt');
savevector3(file,obs,vec.res,vec.res25,vec.res50,vec.res75);
% Save vrig
file = fullfile(savedir,'rig_vector.txt');
savevector(file,obs,vec.rig);
% Save vmec
file = fullfile(savedir,'mec_vector.txt');
savevector3(file,obs,vec.mec,vec.mec25,vec.mec50,vec.mec75);
% Save vkin
file = fullfile(savedir,'kin_vector.txt');
savevector(file,obs,vec.kin);
end

function savevector(file,obs,v)
ve = v(1:3:end);
vn = v(2:3:end);
vu = v(3:3:end);
fid = fopen(file,'wt');
fprintf(fid,'#   site      lon     lat        hei     ve     vn     vu\n');
for nob = 1:obs(1).nobs
  fprintf(fid,'%8s %8.3f %7.3f %10.3f %6.2f %6.2f %6.2f\n',...
      obs(1).name{nob},obs(1).alon(nob),obs(1).alat(nob),obs(1).ahig(nob),...
      ve(nob),vn(nob),vu(nob));
end
fclose(fid);
end

function savevector3(file,obs,v,v1,v2,v3)
ve = v(1:3:end);
vn = v(2:3:end);
vu = v(3:3:end);
ve1 = v1(1:3:end);
vn1 = v1(2:3:end);
vu1 = v1(3:3:end);
ve2 = v2(1:3:end);
vn2 = v2(2:3:end);
vu2 = v2(3:3:end);
ve3 = v3(1:3:end);
vn3 = v3(2:3:end);
vu3 = v3(3:3:end);
fid = fopen(file,'wt');
fprintf(fid,'#   site      lon     lat        hei     ve     vn     vu   ve25   vn25   vu25   ve50   vn50   vu50   ve75   vn75   vu75\n');
for nob = 1:obs(1).nobs
  fprintf(fid,'%8s %8.3f %7.3f %10.3f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n',...
      obs(1).name{nob},obs(1).alon(nob),obs(1).alat(nob),obs(1).ahig(nob),...
      ve(nob),vn(nob),vu(nob),ve1(nob),vn1(nob),vu1(nob),ve2(nob),vn2(nob),vu2(nob),ve3(nob),vn3(nob),vu3(nob));
end
fclose(fid);
end

%% Save coupling and backslip
function SaveBackslip(folder,blk,tcha,s,rep)
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
        rel_st    = s.rel_mec(   mm3m     :mm3m+  nf-1);
        rel_dp    = s.rel_mec(   mm3m+  nf:mm3m+2*nf-1);
        rel_ts    = s.rel_mec(   mm3m+2*nf:mm3m+3*nf-1);
        bslip_st  = s.sdr_mec(   mm3m     :mm3m+  nf-1);
        bslip_dp  = s.sdr_mec(   mm3m+  nf:mm3m+2*nf-1);
        bslip_ts  = s.sdr_mec(   mm3m+2*nf:mm3m+3*nf-1);
        bslip_st25= s.sdr_mec25( mm3m     :mm3m+  nf-1);
        bslip_dp25= s.sdr_mec25( mm3m+  nf:mm3m+2*nf-1);
        bslip_ts25= s.sdr_mec25( mm3m+2*nf:mm3m+3*nf-1);
        bslip_st50= s.sdr_mec50( mm3m     :mm3m+  nf-1);
        bslip_dp50= s.sdr_mec50( mm3m+  nf:mm3m+2*nf-1);
        bslip_ts50= s.sdr_mec50( mm3m+2*nf:mm3m+3*nf-1);
        bslip_st75= s.sdr_mec75( mm3m     :mm3m+  nf-1);
        bslip_dp75= s.sdr_mec75( mm3m+  nf:mm3m+2*nf-1);
        bslip_ts75= s.sdr_mec75( mm3m+2*nf:mm3m+3*nf-1);
        bslipl_st = s.sdr_mecl(mm3m     :mm3m+  nf-1);
        bslipl_dp = s.sdr_mecl(mm3m+  nf:mm3m+2*nf-1);
        bslipl_ts = s.sdr_mecl(mm3m+2*nf:mm3m+3*nf-1);
        id_lock   = s.idl(mm3m     :mm3m+  nf-1);
        id_lock25 = s.idl25(mm3m     :mm3m+  nf-1);
        id_lock50 = s.idl50(mm3m     :mm3m+  nf-1);
        id_lock75 = s.idl75(mm3m     :mm3m+  nf-1);
        slip_st   = s.sr_mec( mm3m     :mm3m+  nf-1);
        slip_dp   = s.sr_mec( mm3m+  nf:mm3m+2*nf-1);
        slip_ts   = s.sr_mec( mm3m+2*nf:mm3m+3*nf-1);
        slipl_st  = s.sr_mecl(mm3m     :mm3m+  nf-1);
        slipl_dp  = s.sr_mecl(mm3m+  nf:mm3m+2*nf-1);
        slipl_ts  = s.sr_mecl(mm3m+2*nf:mm3m+3*nf-1);
        pc = tcha.aveaid(mm1m:mm1m+nf-1,:,rep);
        outdata = [fltnum',...
            blk(1).bound(nb1,nb2).blon,...
            blk(1).bound(nb1,nb2).blat,...
            blk(1).bound(nb1,nb2).bdep,...
            clon, clat, cdep,...
            rel_st, rel_dp, rel_ts, pc,...
            id_lock, bslip_st, bslip_dp, bslip_ts,...
            id_lock25, bslip_st25, bslip_dp25, bslip_ts25,...
            id_lock50, bslip_st50, bslip_dp50, bslip_ts50,...
            id_lock75, bslip_st75, bslip_dp75, bslip_ts75];
        fprintf(fmec,'#    1        2        3        4       5       6       7       8       9      10       11      12      13      14      15      16   17   18      19      20      21     22       23       24       25     26       27       28       29     30       31       32       33\n');
        fprintf(fmec,'#  tri     lon1     lon2     lon3    lat1    lat2    lat3    dep1    dep2    dep3     clon    clat    cdep  rel_st  rel_dp  rel_ts  P_l id_l  sdr_st  sdr_dp  sdr_ts id_l25 sdr_st25 sdr_dp25 sdr_ts25 id_l50 sdr_st50 sdr_dp50 sdr_ts50 id_l75 sdr_st75 sdr_dp75 sdr_ts75\n');
        fprintf(fmec,'%6i %8.3f %8.3f %8.3f %7.3f %7.3f %7.3f %7.2f %7.2f %7.2f %8.3f %7.3f %7.2f %7.2f %7.2f %7.2f %4.2f %4i %7.2f %7.2f %7.2f %6i %8.2f %8.2f %8.2f %6i %8.2f %8.2f %8.2f %6i %8.2f %8.2f %8.2f\n',outdata');
        fclose(fmec);
        mm1m = mm1m +   nf;
        mm3m = mm3m + 3*nf;
      else
        file = fullfile(savedir,['trikin_',num2str(nb1),'_',num2str(nb2),'.txt']);
        fkin = fopen(file,'wt');
        rel_st  = s.rel_kin(mm3k     :mm3k+  nf-1);
        rel_dp  = s.rel_kin(mm3k+  nf:mm3k+2*nf-1);
        rel_ts  = s.rel_kin(mm3k+2*nf:mm3k+3*nf-1);
        sdr_st  = s.sdr_kin(mm3k     :mm3k+  nf-1);
        sdr_dp  = s.sdr_kin(mm3k+  nf:mm3k+2*nf-1);
        sdr_ts  = s.sdr_kin(mm3k+2*nf:mm3k+3*nf-1);
        sr_st   =  s.sr_kin(mm3k     :mm3k+  nf-1);
        sr_dp   =  s.sr_kin(mm3k+  nf:mm3k+2*nf-1);
        sr_ts   =  s.sr_kin(mm3k+2*nf:mm3k+3*nf-1);
        cpmean  = tcha.aveflt(mm1k:mm1k+nf-1,:,rep);
        outdata = [fltnum',...
            blk(1).bound(nb1,nb2).blon,...
            blk(1).bound(nb1,nb2).blat,...
            blk(1).bound(nb1,nb2).bdep,...
            clon, clat, cdep,...
            rel_st, rel_dp, rel_ts,...
            cpmean, sdr_st, sdr_dp, sdr_ts];
        fprintf(fkin,'#    1        2        3        4       5       6       7       8       9      10       11      12      13      14      15      16    14      15      16      17\n');
        fprintf(fkin,'#  tri     lon1     lon2     lon3    lat1    lat2    lat3    dep1    dep2    dep3     clon    clat    cdep  rel_st  rel_dp  rel_ts    CC  sdr_st  sdr_dp  sdr_ts\n');
        fprintf(fkin,'%6i %8.3f %8.3f %8.3f %7.3f %7.3f %7.3f %7.2f %7.2f %7.2f %8.3f %7.3f %7.2f %7.2f %7.2f %7.2f %5.2f %7.2f %7.2f %7.2f\n',outdata');
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
function SaveUDlines(folder,blk)
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
function [s,vec] = CalcOptimumValue(prm,obs,tcha,G,d,rep)
mpmean = tcha.avepol(:,:,rep);
mcmean = tcha.aveflt(:,:,rep);
mamean = tcha.aveasp(:,:,rep);
mimean = tcha.aveine(:,:,rep);
idmean = tcha.aveaid(:,:,rep);

Hu = Heaviside(G(1).zulim - G(1).zc);
Hd = Heaviside(G(1).zdlim - G(1).zc);
Hlim = Hd - Hu;

idl1 = ((Heaviside(G(1).zd*mamean-G(1).zc) - Heaviside(G(1).zu*mamean-G(1).zc)) > 0) .* Hlim;
idl = logical(d(1).maid *  idl1);
idc = logical(d(1).maid * ~idl1);
% -----
% load('depid_5km.mat')  % depth threthould (5km)
% load('depid_6km.mat')  % depth threthould (6km)
% load('depid_8km.mat')  % depth threthould (8km)
% load('depid_10km.mat')  % depth threthould (10km)
% idmean = idmean .* cdep_iddep;  % depth threthould
% load('id_muroto_lock.mat')  % Given locking
% load('lock_id_1968eq.mat')  % 1968 asperity
% idmean(c_idin) = 1;
% load('Result_red/Test_27/Fukushima_30_50km/fukushima_30_50km.mat')   % Given locking
% load('Result_red/Test_27/Miyagi_20_40km/miyagi_20_40km.mat')   % Given locking
% load('Result_red/Test_27/NTohoku_0_10km/ntohoku_0_10km.mat')   % Given locking
% idmean(cdep_iddep) = 1;
% -----
idl25 = logical(d(1).maid *  (idmean>=0.25));
idl50 = logical(d(1).maid *  (idmean>=0.50));
idl75 = logical(d(1).maid *  (idmean>=0.75));
idc25 = logical(d(1).maid * ~(idmean>=0.25));
idc50 = logical(d(1).maid * ~(idmean>=0.50));
idc75 = logical(d(1).maid * ~(idmean>=0.75));

% Calculate back-slip on locked and creeping patches.
rel_mec    = (G(1).tb_mec * mpmean) .* d(1).cfinv_mec;
rel_kin    = (G(1).tb_kin * mpmean) .* d(1).cfinv_kin;
bslip      = rel_mec .* idl;
bslip25    = rel_mec .* idl25;
bslip50    = rel_mec .* idl50;
bslip75    = rel_mec .* idl75;
bslipl     = bslip;
bslipl25   = bslip25;
bslipl50   = bslip50;
bslipl75   = bslip75;

bslip(idc) = -G(1).s(idc,idc) \ (G(1).s(idc,idl) * bslip(idl));
bslip25(idc25) = -G(1).s(idc25,idc25) \ (G(1).s(idc25,idl25) * bslip25(idl25));
bslip50(idc50) = -G(1).s(idc50,idc50) \ (G(1).s(idc50,idl50) * bslip50(idl50));
bslip75(idc75) = -G(1).s(idc75,idc75) \ (G(1).s(idc75,idl75) * bslip75(idl75));
s.idl      = idl;
s.idl25    = idl25;
s.idl50    = idl50;
s.idl75    = idl75;
s.idc      = idc;
s.idc25    = idc25;
s.idc50    = idc50;
s.idc75    = idc75;
s.rel_mec  = rel_mec;
s.rel_kin  = rel_kin;
s.sdr_mec  = bslip;
s.sdr_mec25  = bslip25;
s.sdr_mec50  = bslip50;
s.sdr_mec75  = bslip75;
s.sdr_mecl   = bslipl;
s.sdr_mecl25 = bslipl25;
s.sdr_mecl50 = bslipl50;
s.sdr_mecl75 = bslipl75;
s.sdr_kin    = (s.rel_kin .* (d(1).mcid * mcmean));
s.sr_mec     = -(s.rel_mec - s.sdr_mec );
s.sr_mec25   = -(s.rel_mec - s.sdr_mec25 );
s.sr_mec50   = -(s.rel_mec - s.sdr_mec50 );
s.sr_mec75   = -(s.rel_mec - s.sdr_mec75 );
s.sr_mecl    = -(s.rel_mec - s.sdr_mecl);
s.sr_mecl25  = -(s.rel_mec - s.sdr_mecl25);
s.sr_mecl50  = -(s.rel_mec - s.sdr_mecl50);
s.sr_mecl75  = -(s.rel_mec - s.sdr_mecl75);
s.sr_kin     = -(s.rel_kin - s.sdr_kin );

% Calc vectors for mean parameters
vec.rig = G(1).p * mpmean;
vec.kin = G(1).c_kin * s.sdr_kin;
vec.mec = G(1).c_mec * s.sdr_mec;
vec.mec25 = G(1).c_mec * s.sdr_mec25;
vec.mec50 = G(1).c_mec * s.sdr_mec50;
vec.mec75 = G(1).c_mec * s.sdr_mec75;
vec.ine = G(1).i * mimean;
% Zero padding
if prm.gpu ~= 99
  if isempty(vec.rig); vec.rig = zeros(size(d(1).ind),precision,'gpuArray'); end
  if isempty(vec.kin); vec.kin = zeros(size(d(1).ind),precision,'gpuArray'); end
  if isempty(vec.mec); vec.mec = zeros(size(d(1).ind),precision,'gpuArray'); end
  if isempty(vec.mec25); vec.mec25 = zeros(size(d(1).ind),precision,'gpuArray'); end
  if isempty(vec.mec50); vec.mec50 = zeros(size(d(1).ind),precision,'gpuArray'); end
  if isempty(vec.mec75); vec.mec75 = zeros(size(d(1).ind),precision,'gpuArray'); end
  if isempty(vec.ine); vec.ine = zeros(size(d(1).ind),precision,'gpuArray'); end
else
  if isempty(vec.rig); vec.rig = zeros(size(d(1).ind)); end
  if isempty(vec.kin); vec.kin = zeros(size(d(1).ind)); end
  if isempty(vec.mec); vec.mec = zeros(size(d(1).ind)); end
  if isempty(vec.mec25); vec.mec25 = zeros(size(d(1).ind)); end
  if isempty(vec.mec50); vec.mec50 = zeros(size(d(1).ind)); end
  if isempty(vec.mec75); vec.mec75 = zeros(size(d(1).ind)); end
if isempty(vec.ine); vec.ine = zeros(size(d(1).ind)); end
end
% Total velocities
vec.sum = vec.rig + vec.kin + vec.mec + vec.ine;
vec.sum25 = vec.rig + vec.kin + vec.mec25 + vec.ine;
vec.sum50 = vec.rig + vec.kin + vec.mec50 + vec.ine;
vec.sum75 = vec.rig + vec.kin + vec.mec75 + vec.ine;
vobs = reshape([obs(1).evec; obs(1).nvec; obs(1).hvec],3*obs(1).nobs,1);
vec.res = vobs - vec.sum;
vec.res25 = vobs - vec.sum25;
vec.res50 = vobs - vec.sum50;
vec.res75 = vobs - vec.sum75;

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
