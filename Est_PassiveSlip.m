%% Est_PassiveSlip.m
function Est_PassiveSlip
% Coded by Ryohei Sasajima final 2013/12/23
% Combined by Hiroshi Kimura     2018/11/12
% Revise by Hiroshi Kimura       2019/2/16
%--- 
prm.input = 'PARAMETER/parameter_forward.txt';
prm.optfile='PARAMETER/opt_bound_par_forward.txt';
%--
% Read Parameters.
[prm]     = ReadParameters(prm);
% Read observation data. 
[obs]     = ReadObs(prm);
% Read blocks.
[blk,obs] = ReadBlockBound(prm,obs);
% Read or generate block interfaces.
[blk]     = ReadBlockInterface(blk,prm);
% Read Euler pole vectors.
[eul,prm] = ReadEulerPoles(blk,prm);
% Read dipping boundaries.
[blk,prm] = ReadDippingBound(blk,prm);
% Read internal strain flag.
[blk,prm] = ReadInternalStrain(blk,obs,prm);
% Read locked patches definition.
[blk]     = ReadLockedPatch(blk,prm);
% Estimate rigid motion and calculate AIC.
[blk,obs] = CalcAIC(blk,obs,eul,prm);
% Generate Green's functions.
[tri]     = GreenTri(blk,obs,prm);
% Arrange Green's functions.
[d,G]     = GreenFunctionMatrix(blk,obs,tri);
% Define locking patches.
[d]       = InitialLockingPatch(blk,tri,d);
% Calculate pasive slip and response to surface.
[cal]     = CalcSlip(blk,tri,prm,obs,eul,d,G);
% Save data.
% SaveData(prm,blk,obs,tri,d,G,cal)

end

%% Read Parameter configuration file
function [prm] = ReadParameters(prm)
% MCMC Inversion for Geodetic 
% Coded    by Takeo Ito 2011/11/08 (ver 1.0)
% Modified by Takeo Ito 2012/10/26 (ver 1.1)
% Modified by Takeo Ito 2015/11/11 (ver 1.2)
% Modified by Takeo Ito 2016/07/06 (ver 1.3)
% Modified by Hiroshi Kimura 2018/04/06 (ver 1.4)
% Modified by Hiroshi Kimura 2019/01/28 (ver 2.0)
%
fid = fopen(prm.input,'r');
fileobs            = fscanf(fid,'%s \n',[1,1]); [~] = fgetl(fid);
dirblock           = fscanf(fid,'%s \n',[1,1]); [~] = fgetl(fid);
dirblock_interface = fscanf(fid,'%s \n',[1,1]); [~] = fgetl(fid);
dirblock_patch     = fscanf(fid,'%s \n',[1,1]); [~] = fgetl(fid);
filepole           = fscanf(fid,'%s \n',[1,1]); [~] = fgetl(fid);
filedipb           = fscanf(fid,'%s \n',[1,1]); [~] = fgetl(fid);
fileinternal       = fscanf(fid,'%s \n',[1,1]); [~] = fgetl(fid);
dirresult          = fscanf(fid,'%s \n',[1,1]); [~] = fgetl(fid);
prm.home_d = pwd;
prm.fileobs            = fullfile(prm.home_d,fileobs);
prm.dirblock           = fullfile(prm.home_d,dirblock);
prm.dirblock_interface = fullfile(prm.home_d,dirblock_interface);
prm.dirblock_patch     = fullfile(prm.home_d,dirblock_patch);
prm.filepole           = fullfile(prm.home_d,filepole);
prm.filedipb           = fullfile(prm.home_d,filedipb);
prm.fileinternal       = fullfile(prm.home_d,fileinternal);
prm.dirresult          = fullfile(prm.home_d,dirresult);
%
prm.gpu = fscanf(fid,'%d \n',[1,1]); [~] = fgetl(fid);
prm.itr = fscanf(fid,'%d \n',[1,1]); [~] = fgetl(fid);
prm.thr = fscanf(fid,'%d \n',[1,1]); [~] = fgetl(fid);
prm.cha = fscanf(fid,'%d \n',[1,1]); [~] = fgetl(fid);
prm.kep = fscanf(fid,'%d \n',[1,1]); [~] = fgetl(fid);
prm.rwd = fscanf(fid,'%f \n',[1,1]);
fclose(fid);
%====================================================
tmp = load(prm.optfile);
prm.num    = size(tmp,1);
prm.optb1  = tmp(:,1);
prm.optb2  = tmp(:,2);
prm.optint = tmp(:,3);
%====================================================
fprintf('==================\nINPUT PARAMETERS\n==================\n') 
fprintf('HOME_D                    : %s \n',prm.home_d) 
fprintf('FileOBS                   : %s \n',prm.fileobs) 
fprintf('DIRBlock                  : %s \n',prm.dirblock)
fprintf('DIRBlock_Interface        : %s \n',prm.dirblock_interface) 
fprintf('File fixed epole          : %s \n',prm.filepole) 
fprintf('File Rigid boundary       : %s \n',prm.filedipb) 
fprintf('File Internal deformation : %s \n',prm.fileinternal) 
fprintf('DIRResult                 : %s \n',prm.dirresult) 
fprintf('GPUdev (CPU:99)           : %i \n',prm.gpu) 
fprintf('ITR(Max_Nitr)             : %i \n',prm.itr) 
fprintf('ITR(Threshold_Nitr)       : %i \n',prm.thr) 
fprintf('CHA(Chain)                : %i \n',prm.cha) 
fprintf('KEP(KEEP)                 : %i \n',prm.kep) 
fprintf('RWD(Walk_dis)             : %4.2f \n',prm.rwd) 
fprintf('==================\n') 
%====================================================
disp('PASS READ_PARAMETERS')
end

%% READ OBSERVATION DATA
function [obs] = ReadObs(prm)
%-------------------
% INPUT format Observations:
% unit is mm/yr
% site_name lon lat EW_comp. NS_comp. UD_comp. ERR_EW ERR_NS ERR_UD 
%-------------------
fid_obs = fopen(prm.fileobs,'r');
n = 0;
while 1
  tline = fgetl(fid_obs);
  if ~ischar(tline); break; end
  str   = strsplit(tline);
  n = n+1;
  obs(1).name(n)   = cellstr(str(1));
  obs(1).alon(n)   = str2double(cellstr(str(2)));   % LON
  obs(1).alat(n)   = str2double(cellstr(str(3)));   % LAT
  obs(1).ahig(n)   = str2double(cellstr(str(4)));   % HIG
  obs(1).evec(n)   = str2double(cellstr(str(5)));   % E-W
  obs(1).nvec(n)   = str2double(cellstr(str(6)));   % N-S
  obs(1).hvec(n)   = str2double(cellstr(str(7)));   % U-D
  obs(1).eerr(n)   = str2double(cellstr(str(8)));   % E-W
  obs(1).nerr(n)   = str2double(cellstr(str(9)));   % N-S
  obs(1).herr(n)   = str2double(cellstr(str(10)));  % U-D
  obs(1).weight(n) = str2double(cellstr(str(11)));  % Weight
  obs(1).axyz(n,:) = conv2ell(obs(1).alat(n),obs(1).alon(n),obs(1).ahig(n));
end
obs(1).nobs   = n;
obs(1).ablk   = zeros(obs(1).nobs,1);
obs(1).latmax = max(obs(1).alat);
obs(1).latmin = min(obs(1).alat);
obs(1).lonmax = max(obs(1).alon);
obs(1).lonmin = min(obs(1).alon);
fprintf('==================\n') 
fprintf('Number of GNSS site %4d \n',n)
fprintf('==================\n') 
end

%% Read block boundary data
function [blk,obs] = ReadBlockBound(prm,obs)
ext = '*.txt';
file = dir([prm.dirblock,'/',ext]);
[nblock,~] = size(file);
blk(1).nblock = nblock;
for nb = 1:blk(1).nblock
  tmp = load(fullfile(prm.dirblock,file(nb).name));
  blk(nb).name = file(nb).name;
  blk(nb).lon  = tmp(:,1);
  blk(nb).lat  = tmp(:,2);
end
fprintf('READ BLOCK FILES : %4i \n',blk(1).nblock)
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    blk(1).bound(nb1,nb2).lat = [];
    blk(1).bound(nb1,nb2).lon = [];
    lca = inpolygon(blk(nb1).lon,blk(nb1).lat,blk(nb2).lon,blk(nb2).lat);
    ca  = find(lca);
    if ~isempty(ca) && sum(lca)~=1
      if and(lca(1),lca(end))
        ca0 = find(lca~=true,1,'last')+1:length(lca)-1;
        ca1 = 1:find(lca~=true,1,'first')-1;
        ca  = [ca0 ca1];
      end
      blk(1).bound(nb1,nb2).lat  = blk(nb1).lat(ca);
      blk(1).bound(nb1,nb2).lon  = blk(nb1).lon(ca);
      blk(1).bound(nb1,nb2).bxyz = conv2ell(blk(1).bound(nb1,nb2).lat,blk(1).bound(nb1,nb2).lon,zeros(size(blk(1).bound(nb1,nb2).lon)));
    end
  end
end
%
for n = 1:blk(1).nblock
  ind = inpolygon(obs(1).alon,obs(1).alat,blk(n).lon,blk(n).lat);
  obs(1).ablk(ind) = n;
  obs(n).nblk = sum(ind);
  obs(n).lat  = obs(1).alat(ind);
  obs(n).lon  = obs(1).alon(ind);
  obs(n).hig  = obs(1).ahig(ind);
  obs(n).eve  = obs(1).evec(ind);
  obs(n).nve  = obs(1).nvec(ind);
  obs(n).hve  = obs(1).hvec(ind);
  obs(n).eer  = obs(1).eerr(ind);
  obs(n).ner  = obs(1).nerr(ind);
  obs(n).her  = obs(1).herr(ind);
  obs(n).wgt  = obs(1).weight(ind);
  obs(n).oxyz = conv2ell(obs(n).lat,obs(n).lon,obs(n).hig);
  obs(n).vne  = reshape([obs(1).evec(ind); obs(1).nvec(ind)],obs(n).nblk.*2,1);
  obs(n).ver  = reshape([obs(1).eerr(ind); obs(1).nerr(ind)],obs(n).nblk.*2,1);
  obs(n).vww  = reshape([obs(1).weight(ind)./(obs(1).eerr(ind).^2); obs(1).weight(ind)./(obs(1).nerr(ind).^2)],obs(n).nblk.*2,1);
end
end

%% Read Plate interface
function [blk] = ReadBlockInterface(blk,prm)
% Coded by Takeo Ito 2016/12/21 (ver 1.0)
%
int_tri       =   50;
dep_limit     = -100;
dep_limit_low =  -20;
dirblk        = prm.dirblock_interface;
blk(1).nb = 0;
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    blk(1).bound(nb1,nb2).type = 1;
    pre_tri_f = fullfile(dirblk,['triB_',num2str(nb1),'_',num2str(nb2),'.txt']); 
    fid = fopen(pre_tri_f,'r');
    if fid >= 0
      fprintf('block interface: %2i  %2i \n',nb1,nb2)
      fprintf('read interface tri boudary file : %s \n',pre_tri_f)
      nf   = 0;
      blon = zeros(1,3);
      blat = zeros(1,3);
      bdep = zeros(1,3);
      while 1
        nf    = nf+1;
        loc_f = fscanf(fid,'%f %f %f \n', [3 3]);
        tline = fgetl(fid); if ~ischar(tline); break; end
        blon(nf,:) = loc_f(1,:);  % lon
        blat(nf,:) = loc_f(2,:);  % lat
        bdep(nf,:) = loc_f(3,:);  % hight
        tline = fgetl(fid); if ~ischar(tline); break; end
      end
      fclose(fid);
      bo_tri_f = fullfile(dirblk,['triBO_',num2str(nb1),'_',num2str(nb2),'.txt']); 
      fid = fopen(bo_tri_f,'r');
      if fid >= 0
        bound_blk = textscan(fid,'%f%f'); fclose(fid);
        bound_blk = cell2mat(bound_blk);
        bslon = mean(blon,2);
        bslat = mean(blat,2);
        id = inpolygon(bslon,bslat,bound_blk(:,1),bound_blk(:,2));
        blon = blon(id,:);
        blat = blat(id,:);
        bdep = bdep(id,:);
      end
      blk(1).bound(nb1,nb2).blon = blon;  % lon
      blk(1).bound(nb1,nb2).blat = blat;  % lat
      blk(1).bound(nb1,nb2).bdep = bdep;  % hight
    else
      sub_f = fullfile(dirblk,['B_',num2str(nb1),'_',num2str(nb2),'.txt']);
      fid   = fopen(sub_f,'r');
      if fid >= 0
        fprintf('block interface: %2i  %2i \n',nb1,nb2)
        fprintf('read interface boudary shape file : %s \n',sub_f)
        dep_blk = textscan(fid,'%f%f%f'); fclose(fid);
        dep_blk = cell2mat(dep_blk);
        f    = scatteredinterpolant(dep_blk(:,1),dep_blk(:,2),dep_blk(:,3));
        bo_f = fullfile(dirblk,['BO_',num2str(nb1),'_',num2str(nb2),'.txt']);
        fid  = fopen(bo_f,'r');
        if fid >= 0
          fprintf('read interface boudary definition file : %s \n',bo_f)
          bound_blk = textscan(fid,'%f%f'); fclose(fid);     
          bound_blk = cell2mat(bound_blk);
        else
          idb       = boundary(dep_blk(:,1),dep_blk(:,2));
          bound_blk = dep_blk(idb,:);
        end
        inb = intersect(find(prm.optb1==nb1),find(prm.optb2==nb2));
        if isempty(inb)
          int_bo = int_tri;
        else
          int_bo = prm.optint(inb);
        end
        [p,bstri] = mesh2d_uni(bound_blk,int_bo,bound_blk);
        bslon = p(:,1);
        bslat = p(:,2);
        bsdep = f(bslon,bslat);
       else
        bstri = [];
        leng  = length(blk(1).bound(nb1,nb2).lon);
        if leng~=0
          bslon = [blk(1).bound(nb1,nb2).lon; blk(1).bound(nb1,nb2).lon(1); (blk(1).bound(nb1,nb2).lon(1:leng-1)+blk(1).bound(nb1,nb2).lon(2:leng))./2;blk(1).bound(nb1,nb2).lon(leng)];
          bslat = [blk(1).bound(nb1,nb2).lat; blk(1).bound(nb1,nb2).lat(1); (blk(1).bound(nb1,nb2).lat(1:leng-1)+blk(1).bound(nb1,nb2).lat(2:leng))./2;blk(1).bound(nb1,nb2).lat(leng)];
          bsdep = [zeros(leng,1)            ; dep_limit_low.*ones(leng+1,1)];
          bstri(1:leng-1     ,1:3) = [1     :leng-1;      2:leng    ; leng+2:2*leng]';
          bstri(leng:2*leng-1,1:3) = [leng+1:2*leng; leng+2:2*leng+1;      1:  leng]';
          fprintf('block interface: %2i  %2i auto set %4i \n',nb1,nb2,(leng-1)*2+1)
          blk(1).bound(nb1,nb2).type = 5;
        else
          blk(1).bound(nb1,nb2).blon = [];
          blk(1).bound(nb1,nb2).blat = [];
          blk(1).bound(nb1,nb2).bdep = [];          
        end
      end
      if ~isempty(bstri)
        blk(1).bound(nb1,nb2).blon = [bslon(bstri(:,1)), bslon(bstri(:,2)), bslon(bstri(:,3))];
        blk(1).bound(nb1,nb2).blat = [bslat(bstri(:,1)), bslat(bstri(:,2)), bslat(bstri(:,3))];
        blk(1).bound(nb1,nb2).bdep = [bsdep(bstri(:,1)), bsdep(bstri(:,2)), bsdep(bstri(:,3))];
%
        out_tri_f = fullfile(prm.dirblock,['triB_',num2str(nb1),'_',num2str(nb2),'.out']);
        nlen      = length(blk(1).bound(nb1,nb2).blat(:,1));
        fid_out   = fopen(out_tri_f,'w+');
        fprintf(fid_out,'%10.5f %9.5f %9.3f \n%10.5f %9.5f %9.3f \n%10.5f %9.5f %9.3f \n%10.5f %9.5f %9.3f \n> \n',...
        reshape([blk(1).bound(nb1,nb2).blon(:,1), blk(1).bound(nb1,nb2).blat(:,1), blk(1).bound(nb1,nb2).bdep(:,1),...
                 blk(1).bound(nb1,nb2).blon(:,2), blk(1).bound(nb1,nb2).blat(:,2), blk(1).bound(nb1,nb2).bdep(:,2),...
                 blk(1).bound(nb1,nb2).blon(:,3), blk(1).bound(nb1,nb2).blat(:,3), blk(1).bound(nb1,nb2).bdep(:,3),...
                 blk(1).bound(nb1,nb2).blon(:,1), blk(1).bound(nb1,nb2).blat(:,1), blk(1).bound(nb1,nb2).bdep(:,1)]',4*nlen,3));
        fclose(fid_out);
      end
    end
    blk(1).nb = blk(1).nb+size(blk(1).bound(nb1,nb2).blon,1);
  end
end
end

%% Read locked patches
function [blk] = ReadLockedPatch(blk,prm)
%    Test version coded by H. Kimura 2019/1/29
% Revised version coded by H. Kimura 2019/2/5

for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    blk(1).bound(nb1,nb2).patchid = 0;
    patchfile = fullfile(prm.dirblock_patch,['patchb_',num2str(nb1),'_',num2str(nb2),'.txt']);
    fid       = fopen(patchfile,'r');
    if fid >= 0
      blk(1).bound(nb1,nb2).patchid = 1;
      np    = 0;
      n     = 0;
      while 1
        tline = fgetl(fid);
        if ~ischar(tline) ; break; end
        if tline(1) ~= '>'
          n   = n+1;
          tmp = strsplit(strtrim(tline));
          blk(1).bound(nb1,nb2).patch(np+1).lon(n) = str2double(cellstr(tmp(1)));
          blk(1).bound(nb1,nb2).patch(np+1).lat(n) = str2double(cellstr(tmp(2)));
        else
          np = np+1;
          n  = 0;
          continue;
        end
      end
    end
  end
end

fprintf('=== Read Locked Patches=== \n');
end

%% Save data
function SaveData(prm,blk,obs,tri,d,G,cal)

logfile=fullfile(prm.dirresult,'log.txt');
log_fid=fopen(logfile,'a');
ext=fullfile(prm.dirresult,'Test_*');
file=dir(ext);
if size(file,1)~=0
  d_no=zeros(size(file));
  for ii=1:size(file,1)
    namsplit=strsplit(file(ii).name,'_');
    d_no(ii)=str2double(char(namsplit(2)));
  end
  b=sort(d_no);
  next_no=b(end)+1;
  a_dir=fullfile(prm.dirresult,['Test_',num2str(next_no,'%02i')]);
else
  a_dir=fullfile(prm.dirresult,'Test_01');
end
f_fir=fullfile(a_dir,'figure');
mkdir(a_dir);
mkdir(f_fir);
% 
fprintf('write output file: %s \n',prm.dirresult)
%
save(fullfile(a_dir,'blk.mat'),'blk','-v7.3')
save(fullfile(a_dir,'tri.mat'),'tri','-v7.3')
save(fullfile(a_dir,'prm.mat'),'prm','-v7.3')
save(fullfile(a_dir,'obs.mat'),'obs','-v7.3')
save(fullfile(a_dir,'cal.mat'),'cal','-v7.3')
save(fullfile(a_dir,'grn.mat'),'d','G','-v7.3')
% 
savefig(100,fullfile(f_fir,'back_slip'))
%{
savefig(140,fullfile(f_fir,'vec_rig_ela'))
savefig(130,fullfile(f_fir,'vector'))
savefig(120,fullfile(f_fir,'pole'))
savefig(110,fullfile(f_fir,'std'))
savefig(100,fullfile(f_fir,'coupling'))
%}
fprintf(log_fid,'MODEL= %s\n',prm.dirblock);
fprintf(log_fid,'OBSDATA= %s\n',prm.fileobs);
fclose('all');
movefile(logfile,a_dir);
end

%% Calculate AIC
function [blk,obs] = CalcAIC(blk,obs,eul,prm)
t_sig = 0;
num_b = 0;
blk(1).pole = [];
logfile = fullfile(prm.dirresult,'log.txt');
logfid  = fopen(logfile,'wt');
for n = 1:blk(1).nblock
  sig  = 0;
  evne = [];
  pole = [0; 0; 0];
  obs(n).eev = zeros(obs(n).nblk,1);
  obs(n).env = zeros(obs(n).nblk,1);
  if obs(n).nblk ~= 0
    sig  = 0;
    evne = [0 0];
    if ismember(n,eul.blid)
      po.wx = eul.wx(eul.blid == n);
      po.wy = eul.wy(eul.blid == n);
      po.wz = eul.wz(eul.blid == n);
      [pole,evne,sig] = est_pole_fix(obs(n).oxyz,obs(n).vne,obs(n).ver,po);
      t_sig = t_sig+sig.*2.*obs(n).nblk;
    elseif obs(n).nblk >= 1
      num_b = num_b+1;
      [pole,evne,sig] = est_pole_w(obs(n).oxyz,obs(n).vne,obs(n).vww);
      t_sig = t_sig+sig.*2.*obs(n).nblk;
    end
  end
  blk(n).sig  = sig;
  blk(n).pol  = pole;
  blk(1).pole = [blk(1).pole;blk(n).pol];
  obs(n).eev  = evne(1:2:end);
  obs(n).env  = evne(2:2:end);
  fprintf('BLOCK=%2d NUM_OBS=%2d Sigma^2=%5.2f ',n,obs(n).nblk,sig)
  fprintf(logfid,'BLOCK=%2d NUM_OBS=%2d Sigma^2=%5.2f ',n,obs(n).nblk,sig);
  [latp,lonp,ang] = xyzp2lla(pole(1),pole(2),pole(3));
  fprintf('Lat:%7.2f deg. Lon:%8.2f deg. Ang:%9.2e deg./m.y. \n',latp,lonp,ang);    
  fprintf(logfid,'Lat:%7.2f deg. Lon:%8.2f deg. Ang:%9.2e deg./m.y. \n',latp,lonp,ang);    
end
aic=(obs(1).nobs.*2).*log(t_sig./(obs(1).nobs.*2))+2.*num_b.*3;
caic=aic+2.*num_b.*3.*(num_b.*3+1)./(obs(1).nobs.*2-num_b.*3-1);
fprintf('Sigma^2=%8.3f AIC=%7.3f cAIC=%7.3f K=%2d\n',t_sig./(obs(1).nobs.*2),aic,caic,num_b.*3)
fprintf(logfid,'Sigma^2=%8.3f AIC=%7.3f cAIC=%7.3f K=%2d\n',t_sig./(obs(1).nobs.*2),aic,caic,num_b.*3);
fclose(logfid);
%
end

%% Estimate Euler pole and calculate velocities
function [pl,evne,sigma] = est_pole_w(oxyz,vne,w)
[nobs,~] = size(oxyz);
r        = zeros(nobs.*2,3);
%r(:,1) = -oxyz(:,2).*pvec(3) + pvec(2).*oxyz(:,3);
%r(:,2) = -oxyz(:,3).*pvec(1) + pvec(3).*oxyz(:,1);
%r(:,3) = -oxyz(:,1).*pvec(2) + pvec(1).*oxyz(:,2);
for n = 1:nobs
  r(2.*n-1,1) = -oxyz(n,7).*oxyz(n,3);
  r(2.*n-1,2) = -oxyz(n,5).*oxyz(n,3);
  r(2.*n-1,3) =  oxyz(n,5).*oxyz(n,2)+oxyz(n,7).*oxyz(n,1);
  r(2.*n,1)   =  oxyz(n,4).*oxyz(n,5).*oxyz(n,3)+oxyz(n,6).*oxyz(n,2);
  r(2.*n,2)   = -oxyz(n,4).*oxyz(n,7).*oxyz(n,3)-oxyz(n,6).*oxyz(n,1);
  r(2.*n,3)   =  oxyz(n,4).*oxyz(n,7).*oxyz(n,2)-oxyz(n,4).*oxyz(n,5).*oxyz(n,1);
end
[pl,~,sigma] = lscov(r,vne,w);
evne         = r*pl;
end

%% Calculate velocities from fixed Euler pole
function [pl,evne,sigma] = est_pole_fix(oxyz,vne,w,po)
[nobs,~] = size(oxyz);
r        = zeros(nobs.*2,3);
%r(:,1) = -oxyz(:,2).*pvec(3) + pvec(2).*oxyz(:,3);
%r(:,2) = -oxyz(:,3).*pvec(1) + pvec(3).*oxyz(:,1);
%r(:,3) = -oxyz(:,1).*pvec(2) + pvec(1).*oxyz(:,2);
for n = 1:nobs
  r(2.*n-1,1) = -oxyz(n,7).*oxyz(n,3);
  r(2.*n-1,2) = -oxyz(n,5).*oxyz(n,3);
  r(2.*n-1,3) =  oxyz(n,5).*oxyz(n,2)+oxyz(n,7).*oxyz(n,1);
  r(2.*n,1)   =  oxyz(n,4).*oxyz(n,5).*oxyz(n,3)+oxyz(n,6).*oxyz(n,2);
  r(2.*n,2)   = -oxyz(n,4).*oxyz(n,7).*oxyz(n,3)-oxyz(n,6).*oxyz(n,1);
  r(2.*n,3)   =  oxyz(n,4).*oxyz(n,7).*oxyz(n,2)-oxyz(n,4).*oxyz(n,5).*oxyz(n,1);
end
pl    = [po.wx;po.wy;po.wz];
evne  = r*pl;
sigma = (1/(2*nobs))*sum(((evne-vne)./w).^2);
end

%% Read Euler Pole file
function [eul,prm] = ReadEulerPoles(blk,prm)
% Fix euler poles at the block which has no observation site.
% blid  : block id that includes fix pole
% omega : unit is deg/myr
eul.id   = false;
eul.flag = 0;
eul.blid = 0;
eul.fixw = zeros(3*blk(1).nblock,1);
if exist(prm.filepole,'file')~=2; return; end  % No file.
% 
eul.fixflag=1;
fid    = fopen(prm.filepole,'r');
tmp    = fscanf(fid,'%d %d %f %f %f\n',[5,inf]);
eul.id = zeros(1,blk(1).nblock);
fixrot = zeros(blk(1).nblock,3);
eul.flag  = tmp(1,:);
eul.blid  = tmp(2,:);
eul.lat   = tmp(3,:);
eul.lon   = tmp(4,:);
eul.omega = tmp(5,:);
eul.flag  = logical(eul.flag);   % use or not
eul.blid  = eul.blid(eul.flag) ;
eul.lat   = eul.lat(eul.flag)  ;
eul.lon   = eul.lon(eul.flag)  ;
eul.omega = eul.omega(eul.flag);
eul.lat   = deg2rad(eul.lat)        ;
eul.lon   = deg2rad(eul.lon)        ;
eul.omega = deg2rad(eul.omega.*1e-6);
eul.wx    = eul.omega .* cos(eul.lat) .* cos(eul.lon);
eul.wy    = eul.omega .* cos(eul.lat) .* sin(eul.lon);
eul.wz    = eul.omega .* sin(eul.lat)                ;
eul.id(eul.blid)   = true;
fixrot(eul.blid,1) = eul.wx;
fixrot(eul.blid,2) = eul.wy;
fixrot(eul.blid,3) = eul.wz;
eul.id   = logical(reshape(repmat(eul.id,3,1),3*blk(1).nblock,1));
eul.fixw = reshape(fixrot',3*blk(1).nblock,1);
% 
prm.aprioripole = tmp';
% 
end

%% Read dipping block boundary
function [blk,prm] = ReadDippingBound(blk,prm)
blk(1).dipbo = zeros(1,3);
if exist(prm.filedipb,'file') ~= 2; return; end  % File not exist
fid = fopen(prm.filedipb,'r');
tmp = fscanf(fid,'%d %d %d\n',[3 Inf]);
blk(1).dipbo = tmp';
end

%% Read internal strain file
function [blk,prm] = ReadInternalStrain(blk,obs,prm)
% Coded by Hiroshi Kimura 2018/05/01 (ver 1.0)
% Revised by Hiroshi Kimura 2018/05/08 (ver 1.1)
% Revised by Hiroshi Kimura 2019/01/31 (ver 1.2)
%------------------------- Parameter file format --------------------------
% Block_Num. Flag1
% Flag1 : If uniform internal deformation is calculated, 1, else, 0
% If there is no description about a certain block, internal deformation is
% calculated for the block. If you 'do not' want to calculate internal
% deformation for a certain block, describe the block number and set FLAG1
% to '0'.
% -------------------------------------------------------------------------
blk(1).internal = zeros(1,5);
if exist(prm.fileinternal,'file') ~= 2; return; end
fid = fopen(prm.fileinternal,'r');
tmp = fscanf(fid,'%d %d\n',[2 inf]);
blk(1).internal = tmp';
idinter = zeros(1,blk(1).nblock);
for nb = 1:blk(1).nblock
  [id,flag] = find(blk(1).internal(:,1) == nb);
  blk(nb).flaginter = 0;
  blk(nb).latinter  = mean(obs(nb).lat);
  blk(nb).loninter  = mean(obs(nb).lon);
  if (flag == 1 && blk(1).internal(id,2) == 1) || flag == 0
    blk(nb).flaginter = 1;
  end
  [blk(nb).xinter,blk(nb).yinter] = PLTXY(obs(nb).lat,obs(nb).lon,blk(nb).latinter,blk(nb).loninter);
  idinter(nb) = blk(nb).flaginter;
end
blk(1).idinter = reshape(repmat(idinter,3,1),3*blk(1).nblock,1);
end

%% CalcTriDisps.m
function [U] = CalcTriDisps(sx, sy, sz, x, y, z, pr, ss, ts, ds)
% CalcTriDisps.m
%
% Calculates displacements due to slip on a triangular dislocation in an
% elastic half space utilizing the Comninou and Dunders (1975) expressions
% for the displacements due to an angular dislocation in an elastic half
% space.
%
% Arguments
%  sx : x-coordinates of observation points
%  sy : y-coordinates of observation points
%  sz : z-coordinates of observation points
%  x  : x-coordinates of triangle vertices.
%  y  : y-coordinates of triangle vertices.
%  z  : z-coordinates of triangle vertices.
%  pr : Poisson's ratio
%  ss : strike slip displacement
%  ts : tensile slip displacement
%  ds : dip slip displacement
%
% Returns
%  U  : structure containing the displacements (U.x, U.y, U.z)
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


% Calculate the slip vector in XYZ coordinates
normVec                      = cross([x(2);y(2);z(2)]-[x(1);y(1);z(1)], [x(3);y(3);z(3)]-[x(1);y(1);z(1)]);
normVec                      = normVec./norm(normVec);

[n1,n2]                      = size(sz);
sz1                          = sz(:);

if (normVec(3) < 0) % Enforce clockwise circulation
    normVec                   = -normVec;
    [x(2),x(3)]               = swap(x(2), x(3));
    [y(2),y(3)]               = swap(y(2), y(3));
    [z(2),z(3)]               = swap(z(2), z(3));
end
strikeVec                    = [-sin(atan2(normVec(2),normVec(1))) cos(atan2(normVec(2),normVec(1))) 0];
dipVec                       = cross(normVec, strikeVec);
slipComp                     = [ss ds ts];
slipVec                      = [strikeVec(:) dipVec(:) normVec(:)] * slipComp(:);

% Solution vectors
U.x                          = zeros(size(sx));
U.y                          = zeros(size(sx));
U.z                          = zeros(size(sx));

% Add a copy of the first vertex to the vertex list for indexing
x(4)                         = x(1);
y(4)                         = y(1);
z(4)                         = z(1);

for iTri = 1:3
    % Calculate strike and dip of current leg
    strike                   = 180/pi*(atan2(y(iTri+1)-y(iTri), x(iTri+1)-x(iTri)));
    segMapLength             = sqrt((x(iTri)-x(iTri+1))^2 + (y(iTri)-y(iTri+1))^2);
    [rx,ry]                  = RotateXyVec(x(iTri+1)-x(iTri), y(iTri+1)-y(iTri), -strike);
    dip                      = 180/pi*(atan2(z(iTri+1)-z(iTri), rx));
    
    if dip >= 0
        beta                  = pi/180*(90-dip);
        if beta > pi/2
            beta              = pi/2-beta;
        end
    else
        beta                  = -pi/180*(90+dip);
        if beta < -pi/2
            beta = pi/2-abs(beta);
        end
    end
    
    ssVec                    = [cos(strike/180*pi) sin(strike/180*pi) 0];
    tsVec                    = [-sin(strike/180*pi) cos(strike/180*pi) 0];
    dsVec                    = cross(ssVec, tsVec);
    lss                      = dot(slipVec, ssVec);
    lts                      = dot(slipVec, tsVec);
    lds                      = dot(slipVec, dsVec);
    
    
    if (abs(beta) > 0.000001) && (abs(beta-pi) > 0.000001)
        % First angular dislocation
        [sx1, sy1]                = RotateXyVec(sx-x(iTri), sy-y(iTri), -strike);
        sx1(abs(sx1)<0.0001)=0.0001;%for eliminate NAN by R. Sasajima and T. Ito ,Nagoya. U. in 2012%
        sy1(abs(sy1)<0.0001)=0.0001;
        [ux1, uy1, uz1]           = adv(sx1, sy1, sz1-z(iTri), z(iTri), beta, pr, lss, lts, lds);
        
        % Second angular dislocation
        [sx2, sy2]                = RotateXyVec(sx-x(iTri+1), sy-y(iTri+1), -strike);
        sx2(abs(sx2)<0.0001)=0.0001;%for eliminate NAN by R. Sasajima and T. Ito ,Nagoya. U. in 2012%
        sy2(abs(sy2)<0.0001)=0.0001;
        [ux2, uy2, uz2]           = adv(sx2, sy2, sz1-z(iTri+1), z(iTri+1), beta, pr, lss, lts, lds);
        
        % Rotate vectors to correct for strike
        [uxn, uyn]                = RotateXyVec(ux1-ux2, uy1-uy2, strike);
        uzn                       = uz1-uz2;
        
        % Add the displacements from current leg
        U.x                       = U.x + reshape(uxn,n1,n2);
        U.y                       = U.y + reshape(uyn,n1,n2);
        U.z                       = U.z + reshape(uzn,n1,n2);
    end
end

% Identify indices for stations under current triangle
inPolyIdx                       = find(inpolygon(sx, sy, x, y) == 1);
underIdx = [];
for iIdx = 1 : numel(inPolyIdx)
    d                           = LinePlaneIntersect(x, y, z, sx(inPolyIdx(iIdx)), sy(inPolyIdx(iIdx)), sz(inPolyIdx(iIdx)));
    if d(3)-sz(inPolyIdx(iIdx)) < 0
        underIdx = [underIdx ; inPolyIdx(iIdx)];
    end
end

% Apply static offset to the points that lie underneath the current triangle
U.x(underIdx)                = U.x(underIdx) - slipVec(1);
U.y(underIdx)                = U.y(underIdx) - slipVec(2);
U.z(underIdx)                = U.z(underIdx) - slipVec(3);
end

function d = LinePlaneIntersect(x, y, z, sx, sy, sz)
% Calculate the intersection of a line and a plane using a parametric
% representation of the plane.  This is hardcoded for a vertical line.
numerator                       = [1 1 1 1 ; x(1) x(2) x(3) sx ; y(1) y(2) y(3) sy ; z(1) z(2) z(3) sz];
numerator                       = det(numerator);
denominator                     = [1 1 1 0 ; x(1) x(2) x(3) 0 ; y(1) y(2) y(3) 0 ; z(1) z(2) z(3) -sz];
denominator                     = det(denominator);
if denominator == 0
    denominator                 = eps;
end
t                               = numerator/denominator; % parametric curve parameter
d                               = [sx sy sz]-([sx sy 0]-[sx sy sz])*t;
end

function [xp, yp] = RotateXyVec(x, y, alpha)
% Rotate a vector by an angle alpha
x                             = x(:);
y                             = y(:);
alpha                         = pi/180*alpha;
xp                            = cos(alpha).*x - sin(alpha).*y;
yp                            = sin(alpha).*x + cos(alpha).*y;
end

%% CalcTriStrains.m
function [S] = CalcTriStrains(sx, sy, sz, x, y, z, pr, ss, ts, ds)
% CalcTriStrains.m
%
% Calculates strains due to slip on a triangular dislocation in an
% elastic half space utilizing the symbolically differentiated
% displacement gradient tensor derived from the expressions for
% the displacements due to an angular dislocation in an elastic half
% space (Comninou and Dunders, 1975).
%
% Arguments
%  sx : x-coordinates of observation points
%  sy : y-coordinates of observation points
%  sz : z-coordinates of observation points
%  x  : x-coordinates of triangle vertices.
%  y  : y-coordinates of triangle vertices.
%  z  : z-coordinates of triangle vertices.
%  pr : Poisson's ratio
%  ss : strike slip displacement
%  ts : tensile slip displacement
%  ds : dip slip displacement
%
% Returns
%  S  : structure containing the strains (S.xx, S.yy, S.zz, S.xy, S.xz, S.yz)
%
% This paper should and related code should be cited as:
% Brendan J. Meade, Algorithms for the calculation of exact
% displacements, strains, and stresses for Triangular Dislocation
% Elements in a uniform elastic half space, Computers &
% Geosciences (2007), doi:10.1016/j.cageo.2006.12.003.
%
% Use at your own risk and please let me know of any bugs/errors.
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


% Calculate the slip vector in XYZ coordinates
normVec                      = cross([x(2);y(2);z(2)]-[x(1);y(1);z(1)], [x(3);y(3);z(3)]-[x(1);y(1);z(1)]);
normVec                      = normVec./norm(normVec);
if (normVec(3) < 0) % Enforce clockwise circulation
    normVec                   = -normVec;
    [x(2), x(3)]              = swap(x(2), x(3));
    [y(2), y(3)]              = swap(y(2), y(3));
    [z(2), z(3)]              = swap(z(2), z(3));
end
strikeVec                    = [-sin(atan2(normVec(2),normVec(1))) cos(atan2(normVec(2),normVec(1))) 0];
dipVec                       = cross(normVec, strikeVec);
slipComp                     = [ss ds ts];
slipVec                      = [strikeVec(:) dipVec(:) normVec(:)] * slipComp(:);

% Solution vectors
S.xx                         = zeros(size(sx));
S.yy                         = zeros(size(sx));
S.zz                         = zeros(size(sx));
S.xy                         = zeros(size(sx));
S.xz                         = zeros(size(sx));
S.yz                         = zeros(size(sx));

% Add a copy of the first vertex to the vertex list for indexing
x(4)                         = x(1);
y(4)                         = y(1);
z(4)                         = z(1);

for iTri = 1:3
    % Calculate strike and dip of current leg
    strike                   = 180/pi*(atan2(y(iTri+1)-y(iTri), x(iTri+1)-x(iTri)));
    segMapLength             = sqrt((x(iTri)-x(iTri+1))^2 + (y(iTri)-y(iTri+1))^2);
    [rx, ry]                 = RotateXyVec(x(iTri+1)-x(iTri), y(iTri+1)-y(iTri), -strike);
    dip                      = 180/pi*(atan2(z(iTri+1)-z(iTri), rx));
    
    if dip >= 0
        beta                  = pi/180*(90-dip);
        if beta > pi/2
            beta               = pi/2-beta;
        end
    else
        beta                  = -pi/180*(90+dip);
        if beta < -pi/2
            beta = pi/2-abs(beta);
        end
    end
    ssVec                    = [cos(strike/180*pi) sin(strike/180*pi) 0];
    tsVec                    = [-sin(strike/180*pi) cos(strike/180*pi) 0];
    dsVec                    = cross(ssVec, tsVec);
    lss                      = dot(slipVec, ssVec);
    lts                      = dot(slipVec, tsVec);
    lds                      = dot(slipVec, dsVec);
    
    if (abs(beta) > 0.000001) && (abs(beta-pi) > 0.000001)
        % First angular dislocation
        [sx1,sy1]                 = RotateXyVec(sx-x(iTri), sy-y(iTri), -strike);
        sx1(abs(sx1)<0.0001)=0.0001;%for eliminate NAN by R. Sasajima and T. Ito ,Nagoya. U. in 2012%
        sy1(abs(sy1)<0.0001)=0.0001;
        [a11,a22,a33,a12,a13,a23] = advs(sx1, sy1, sz-z(iTri), z(iTri), beta, pr, lss, lts, lds);
        
        % Second angular dislocation
        [sx2,sy2]                 = RotateXyVec(sx-x(iTri+1), sy-y(iTri+1), -strike);
        sx2(abs(sx2)<0.0001)=0.0001;
        sy2(abs(sy2)<0.0001)=0.0001;
        [b11,b22,b33,b12,b13,b23] = advs(sx2, sy2, sz-z(iTri+1), z(iTri+1), beta, pr, lss, lts, lds);
        
        % Rotate tensors to correct for strike
        bxx                       = a11-b11;
        byy                       = a22-b22;
        bzz                       = a33-b33;
        bxy                       = a12-b12;
        bxz                       = a13-b13;
        byz                       = a23-b23;
        
        g                         = pi/180*strike;
        e11n                      = (cos(g)*bxx-sin(g)*bxy)*cos(g)-(cos(g)*bxy-sin(g)*byy)*sin(g);
        e12n                      = (cos(g)*bxx-sin(g)*bxy)*sin(g)+(cos(g)*bxy-sin(g)*byy)*cos(g);
        e13n                      = cos(g)*bxz-sin(g)*byz;
        e22n                      = (sin(g)*bxx+cos(g)*bxy)*sin(g)+(sin(g)*bxy+cos(g)*byy)*cos(g);
        e23n                      = sin(g)*bxz+cos(g)*byz;
        e33n                      = bzz;
        
        % Add the strains from current leg
        S.xx                      = S.xx + e11n;
        S.yy                      = S.yy + e22n;
        S.zz                      = S.zz + e33n;
        S.xy                      = S.xy + e12n;
        S.xz                      = S.xz + e13n;
        S.yz                      = S.yz + e23n;
    end
end
end

function [a, b] = swap(a, b)
% Swap two values
temp                            = a;
a                               = b;
b                               = temp;
end

%% Make Green's function of halfspace elastic strain for triangular meshes
function [tri] = GreenTri(blk,obs,prm)
% Coded by Hiroshi Kimura 2019/01/28 (ver 1.0)
pr = 0.25;
nd = size(obs(1).alat,2);
%
alat = mean(obs(1).alat(:));
alon = mean(obs(1).alon(:));
[obsx,obsy] = PLTXY(obs(1).alat,obs(1).alon,alat,alon);
obsz = -1e-3.*obs(1).ahig;
%
tri(1).obsdis= [];
tri(1).nb    = 0;
for nb1 = 1:blk(1).nblock
  [blk(nb1).localx,blk(nb1).localy] = PLTXY(blk(nb1).lat,blk(nb1).lon,alat,alon);
  for nb2 = nb1+1:blk(1).nblock
    nf = size(blk(1).bound(nb1,nb2).blon,1);
    if nf ~= 0
      trimat = fullfile(prm.dirblock_interface,['tri_',num2str(nb1),'_',num2str(nb2),'.mat']);
      fid = fopen(trimat);
      if fid > 0
        load(trimat);
        tri(1).bound(nb1,nb2) = trisave;
      else
        tri(1).bound(nb1,nb2).gustr = zeros(3*nd,nf);
        tri(1).bound(nb1,nb2).gudip = zeros(3*nd,nf);
        tri(1).bound(nb1,nb2).gutns = zeros(3*nd,nf);
        tri(1).bound(nb1,nb2).cf   =  ones(3*nf,1);
        tri(1).bound(nb1,nb2).inv  = zeros(3*nf,1);
        %
        fprintf('==================\n Block %2i : Block %2i \n Number of TRI sub-faults : %4i \n',nb1,nb2,nf)
        %
        triclon = mean(blk(1).bound(nb1,nb2).blon,2);
        triclat = mean(blk(1).bound(nb1,nb2).blat,2);
        tricdep = mean(blk(1).bound(nb1,nb2).bdep,2);
        [tricx,tricy] = PLTXY(triclat,triclon,alat,alon);
        tricz = -tricdep;
        %
        blk(1).bound(nb1,nb2).flag1 = 0;
        for ii = 1:size(blk(1).dipbo,1)
          dippingid = ismember([nb1 nb2],blk(1).dipbo(ii,2:3));
          ispair    = sum(dippingid);
          if ispair == 2
            if max(blk(1).dipbo(ii,2:3)) == blk(1).dipbo(ii,1)
              blk(1).bound(nb1,nb2).flag1 = 1;
            else
              blk(1).bound(nb1,nb2).flag1 = 2;
            end
            break
          end
        end
        %
        for n = 1:nf
          [trix,triy] = PLTXY(blk(1).bound(nb1,nb2).blat(n,:),blk(1).bound(nb1,nb2).blon(n,:),alat,alon);
          triz        = -1.*blk(1).bound(nb1,nb2).bdep(n,:);
          f_loc       = [trix;triy;triz];
          [f,da,nv,st,dp,phi,theta] = EstFaultTri(f_loc);
          nv(3) = -nv(3);
          dp(3) = -dp(3);
          tri(1).bound(nb1,nb2).clat(n)   = mean(blk(1).bound(nb1,nb2).blat(n,:));
          tri(1).bound(nb1,nb2).clon(n)   = mean(blk(1).bound(nb1,nb2).blon(n,:));
          tri(1).bound(nb1,nb2).cdep(n)   = mean(blk(1).bound(nb1,nb2).bdep(n,:)); % up is plus
          tri(1).bound(nb1,nb2).da(n)     = da;
          tri(1).bound(nb1,nb2).nv(n,:)   = nv;
          tri(1).bound(nb1,nb2).st(n,:)   = st;
          tri(1).bound(nb1,nb2).dp(n,:)   = dp;
          tri(1).bound(nb1,nb2).phi(n)    = phi;     % strike from x-axis (counter-clock wise)
          tri(1).bound(nb1,nb2).theta(n)  = theta;   % dip from x-y plane (clock wise)
          tri(1).bound(nb1,nb2).oxyz(n,:) = conv2ell(tri(1).bound(nb1,nb2).clat(n),tri(1).bound(nb1,nb2).clon(n),1e3.*tri(1).bound(nb1,nb2).cdep(n));
          u = CalcTriDisps(obsx,obsy,obsz,trix,triy,triz,pr,1,0,0);
          s = CalcTriStrains(tricx,tricy,tricz,trix,triy,triz,pr,1,0,0);
          tri(1).bound(nb1,nb2).gustr(1:3:3*nd,n) =  u.x;  % E
          tri(1).bound(nb1,nb2).gustr(2:3:3*nd,n) =  u.y;  % N
          tri(1).bound(nb1,nb2).gustr(3:3:3*nd,n) = -u.z;  % U
          tri(1).bound(nb1,nb2).gsstr(1:6:6*nf,n) =  s.xx;  % exx
          tri(1).bound(nb1,nb2).gsstr(2:6:6*nf,n) =  s.xy;  % exy
          tri(1).bound(nb1,nb2).gsstr(3:6:6*nf,n) = -s.xz;  % exz
          tri(1).bound(nb1,nb2).gsstr(4:6:6*nf,n) =  s.yy;  % eyy
          tri(1).bound(nb1,nb2).gsstr(5:6:6*nf,n) = -s.yz;  % eyz
          tri(1).bound(nb1,nb2).gsstr(6:6:6*nf,n) = -s.zz;  % ezz
          u = CalcTriDisps(obsx,obsy,obsz,trix,triy,triz,pr,0,1,0);
          s = CalcTriStrains(tricx,tricy,tricz,trix,triy,triz,pr,0,1,0);
          tri(1).bound(nb1,nb2).gutns(1:3:3*nd,n) =  u.x;  % E
          tri(1).bound(nb1,nb2).gutns(2:3:3*nd,n) =  u.y;  % N
          tri(1).bound(nb1,nb2).gutns(3:3:3*nd,n) = -u.z;  % U 
          tri(1).bound(nb1,nb2).gstns(1:6:6*nf,n) =  s.xx;  % exx
          tri(1).bound(nb1,nb2).gstns(2:6:6*nf,n) =  s.xy;  % exy
          tri(1).bound(nb1,nb2).gstns(3:6:6*nf,n) = -s.xz;  % exz
          tri(1).bound(nb1,nb2).gstns(4:6:6*nf,n) =  s.yy;  % eyy
          tri(1).bound(nb1,nb2).gstns(5:6:6*nf,n) = -s.yz;  % eyz
          tri(1).bound(nb1,nb2).gstns(6:6:6*nf,n) = -s.zz;  % ezz
          u = CalcTriDisps(obsx,obsy,obsz,trix,triy,triz,pr,0,0,1);
          s = CalcTriStrains(tricx,tricy,tricz,trix,triy,triz,pr,0,0,1);
          tri(1).bound(nb1,nb2).gudip(1:3:3*nd,n) =  u.x;  % E
          tri(1).bound(nb1,nb2).gudip(2:3:3*nd,n) =  u.y;  % N
          tri(1).bound(nb1,nb2).gudip(3:3:3*nd,n) = -u.z;  % U
          tri(1).bound(nb1,nb2).gsdip(1:6:6*nf,n) =  s.xx;  % exx
          tri(1).bound(nb1,nb2).gsdip(2:6:6*nf,n) =  s.xy;  % exy
          tri(1).bound(nb1,nb2).gsdip(3:6:6*nf,n) = -s.xz;  % exz
          tri(1).bound(nb1,nb2).gsdip(4:6:6*nf,n) =  s.yy;  % eyy
          tri(1).bound(nb1,nb2).gsdip(5:6:6*nf,n) = -s.yz;  % eyz
          tri(1).bound(nb1,nb2).gsdip(6:6:6*nf,n) = -s.zz;  % ezz
          if mod(n,ceil(nf/3)) == 1
            fprintf('MAKE GREEN at TRI sub-faults : %4i / %4i \n',n,nf)
          end
          [tri]     = CorrectFactor(blk,tri,nb1,nb2,dp,n,nf);
          [blk,tri] = DiscriminateDirection(blk,tri,nb1,nb2,trix,triy,n,nf);
          [tri]     = trans_xyz2strdip(tri,nb1,nb2,n,nf);
        end
        trisave = tri(1).bound(nb1,nb2);
        save(trimat,'trisave','-v7.3');
      end
      tri(1).nb = tri(1).nb+nf;
    else
      tri(1).bound(nb1,nb2).clat = [];
      tri(1).bound(nb1,nb2).clon = [];
      tri(1).bound(nb1,nb2).cdep = [];
    end
  end
end
disp('==================')
disp('PASS GREEN_TRI')
disp('==================')
end

%% Calculate correction factor of (STR, DIP, TNS) unit vectors
function [tri] = CorrectFactor(blk,tri,nb1,nb2,dp,n,nf)
% Coded by H.Kimura 2018/1/31 (test ver.)
switch blk(1).bound(nb1,nb2).flag1
  case {1,2}
    cf = 1/sqrt(dp(1)^2+dp(2)^2);      % 1=sqrt(dp(1)^2+dp(2)^2+dp(3)^2): norm of dp
    if cf==Inf, cf=1; end
    tri(1).bound(nb1,nb2).cf(nf+n) = cf;
  case 0
    return
end
end

%% DISCRIMINATE BOUNDARY TYPE AND SUBFAULT SURFACE DIRECTION
function [blk,tri] = DiscriminateDirection(blk,tri,nb1,nb2,trix,triy,n,nf)
% Coded by H.Kimura 2017/4/28 (test ver.)
% Modified by H.Kimura 2018/2/6
switch blk(1).bound(nb1,nb2).flag1
  case 1
    tri(1).bound(nb1,nb2).inv(     n) = 1;
    tri(1).bound(nb1,nb2).inv(  nf+n) = 1;
    tri(1).bound(nb1,nb2).inv(2*nf+n) = 0;    
  case 2
    tri(1).bound(nb1,nb2).inv(     n) = -1;
    tri(1).bound(nb1,nb2).inv(  nf+n) = -1;
    tri(1).bound(nb1,nb2).inv(2*nf+n) =  0;    
  case 0
    trixc   = mean(trix);
    triyc   = mean(triy);
    [in,on] = inpolygon(trixc,triyc,blk(nb1).localx,blk(nb1).localy);
    if in==1 && on~=1
      tri(1).bound(nb1,nb2).inv(     n) = 1;
      tri(1).bound(nb1,nb2).inv(  nf+n) = 0;
      tri(1).bound(nb1,nb2).inv(2*nf+n) = 1;
    elseif in==1 && on==1
      tric = [trixc triyc 0];
      uv   = [0 0 1];
      nv   = cross(uv,tri(1).bound(nb1,nb2).st(n,:));
      cnv  = tric+nv;
      if inpolygon(cnv(1),cnv(2),blk(nb1).localx,blk(nb1).localy)==1
        tri(1).bound(nb1,nb2).inv(     n) = 1;
        tri(1).bound(nb1,nb2).inv(  nf+n) = 0;
        tri(1).bound(nb1,nb2).inv(2*nf+n) = 1;
      else
        tri(1).bound(nb1,nb2).inv(     n) = -1;
        tri(1).bound(nb1,nb2).inv(  nf+n) =  0;
        tri(1).bound(nb1,nb2).inv(2*nf+n) = -1;
      end
    else
      tri(1).bound(nb1,nb2).inv(     n) = -1;
      tri(1).bound(nb1,nb2).inv(  nf+n) =  0;
      tri(1).bound(nb1,nb2).inv(2*nf+n) = -1;
    end
end
% 
end

%% ESTIMATE FAULT PARAMETERS FOR TRI
function [floc,da,nv,st,dp,phi,theta] = EstFaultTri(loc_f)
% Coded by Takeo Ito 2015/11/11 (ver 1.0)
% Modified by Hiroshi Kimura 2019/01/29
[da]       = AreaTri(loc_f(1,:),loc_f(2,:),loc_f(3,:));
floc       = mean(loc_f,2)';
[nv,st,dp,phi,theta] = Tri_StrDip(loc_f(1,:),loc_f(2,:),loc_f(3,:));
end

%% ESTIMATE AREA AT SUB-FAULT FOR TRI
function [da] = AreaTri(x,y,z)
% CALC. AREA IN THREE DIMENSION USING HERON'S FOMULA
% Coded by Takeo Ito 2006/03/04 (ver 1.0)
leng(1) = sqrt((x(1)-x(2)).^2 + (y(1)-y(2)).^2 + (z(1)-z(2)).^2);
leng(2) = sqrt((x(2)-x(3)).^2 + (y(2)-y(3)).^2 + (z(2)-z(3)).^2);
leng(3) = sqrt((x(3)-x(1)).^2 + (y(3)-y(1)).^2 + (z(3)-z(1)).^2);
s1 = (leng(1)+leng(2)+leng(3))./2;
da = sqrt(s1*(s1-leng(1))*(s1-leng(2))*(s1-leng(3)));
end

%% Estimate strike and dip for tri fault
function [nv,st,dp,phi,theta] = Tri_StrDip(x,y,z)
%==========
% Calculate strike and dip on fault.
% Coded by Hiroshi Kimura (2019/01/29)
% Up is minus.
% phi   : Strike angle of from x-axis (counter-clock wise).
% theta : Dip angle from x-y plane.
%==========
nv = cross([x(2);y(2);z(2)]-[x(1);y(1);z(1)], [x(3);y(3);z(3)]-[x(1);y(1);z(1)]);
nv = nv./norm(nv);
if (nv(3) < 0); nv = -nv; end % Enforce clockwise circulation
st = [-sin(atan2(nv(2),nv(1))) cos(atan2(nv(2),nv(1))) 0];
dp = cross(nv,st);
phi   = atan2(st(2),st(1));
theta = asin(dp(3));
end

%%
function [tri] = trans_xyz2strdip(tri,nb1,nb2,n,nf)
% This function transforms a strain tensor from xyz to fault strike-dip.
% Strain of strike direction on the fault corresponds to ezx',
% strain of dip direction on the fault corresponds to ezy', and
% strain of tensile direction on the fault corresponds to ezz'.

c1=cos(tri(1).bound(nb1,nb2).phi(n));
s1=sin(tri(1).bound(nb1,nb2).phi(n));
c2=cos(-tri(1).bound(nb1,nb2).theta(n));
s2=sin(-tri(1).bound(nb1,nb2).theta(n));

[sst,sdp,sts] = calctrans(tri(1).bound(nb1,nb2).gsstr(:,n),c1,s1,c2,s2);  % response to strike slip
tri(1).bound(nb1,nb2).gsstrT(1:3:3*nf,n)=sst;  % str
tri(1).bound(nb1,nb2).gsstrT(2:3:3*nf,n)=sdp;  % dip
tri(1).bound(nb1,nb2).gsstrT(3:3:3*nf,n)=sts;  % tns
[sst,sdp,sts] = calctrans(tri(1).bound(nb1,nb2).gstns(:,n),c1,s1,c2,s2);  % response to tensile slip
tri(1).bound(nb1,nb2).gstnsT(1:3:3*nf,n)=sst;  % str
tri(1).bound(nb1,nb2).gstnsT(2:3:3*nf,n)=sdp;  % dip
tri(1).bound(nb1,nb2).gstnsT(3:3:3*nf,n)=sts;  % tns
[sst,sdp,sts] = calctrans(tri(1).bound(nb1,nb2).gsdip(:,n),c1,s1,c2,s2);  % response to dip slip
tri(1).bound(nb1,nb2).gsdipT(1:3:3*nf,n)=sst;  % str
tri(1).bound(nb1,nb2).gsdipT(2:3:3*nf,n)=sdp;  % dip
tri(1).bound(nb1,nb2).gsdipT(3:3:3*nf,n)=sts;  % tns

end

%%
function [sst,sdp,sts]=calctrans(sxyz,c1,s1,c2,s2)
% E_stdp = Tx * Tz * E_xyz * Tz' * Tx'
% <->
% |exx' exy' exz'| |1   0  0| | c1 s1 0| |exx exy exz| | c1 s1 0|T |1   0  0|T
% |eyx' eyy' eyz'|=|0  c2 s2|*|-s1 c1 0|*|eyx eyy eyz|*|-s1 c1 0| *|0  c2 s2|
% |ezx' ezy' ezz'| |0 -s2 c2| |  0  0 1| |ezx ezy ezz| |  0  0 1|  |0 -s2 c2|
% where 
% c1 : cos(strike)
% s1 : sin(strike)
% c2 : cos(dip)
% s2 : sin(dip)
% Note strike is the counter clock wise angle from x-axis, dip is angle
% from xy-plane. "T" indicates the transpose of matrix.

c1_2 = c1.^2;
c2_2 = c2.^2;
s1_2 = s1.^2;
s2_2 = s2.^2;
c1s1 = c1.*s1;
c2s2 = c2.*s2;
% ezx' = sst = -s2*( c1*s1*( eyy - exx ) + ( c1^2 - c2^2 )*exy ) + c2( c1*exz + s1*eyz );
% ezy' = sdp = c2*s2*( ezz - s1^2*exx + 2*c1*s1*exy - c1^2*eyy ) + ( c2^2 - s2^2 )*( c1*eyz - s1*exz );
% ezz' = sts = s2^2*( s1^2*exx -2*c1*s1*exy + c1^2*eyy ) + 2*c2*s2*( s1*exz - c1*eyz ) + c2^2*ezz;
exx = sxyz(1:6:end,:);
exy = sxyz(2:6:end,:);
exz = sxyz(3:6:end,:);
eyy = sxyz(4:6:end,:);
eyz = sxyz(5:6:end,:);
ezz = sxyz(6:6:end,:);
sst = -s2.*( c1s1.*( eyy - exx ) + ( c1_2 - c2_2 ).*exy )...
         + c2.*( c1.*exz + s1.*eyz );
sdp =  c2.*s2.*( ezz - s1_2.*exx + 2.*c1s1.*exy - c1_2.*eyy )...
         + ( c2_2 - s2_2 ).*( c1.*eyz - s1.*exz );
sts = s2_2.*( s1_2.*exx -2*c1s1.*exy + c1_2.*eyy )...
         + 2.*c2s2.*( s1.*exz - c1.*eyz )...
         + c2_2.*ezz;
end

%% Make Matrix
function [d,G] = GreenFunctionMatrix(blk,obs,tri)
% Coded by Takeo Ito 2017/01/02 (ver 1.1)
% Coded by Hiroshi Kimura 2018/05/01 (ver 1.2)
% Modified by Hiroshi Kimura 2019/01/30 (ver 1.3)
% pole unit is mm
nobs = length(obs(1).evec);
tmp.obs(1:3:3*nobs) = obs(1).evec;
tmp.obs(2:3:3*nobs) = obs(1).nvec;
tmp.obs(3:3:3*nobs) = obs(1).hvec;
tmp.err(1:3:3*nobs) = obs(1).eerr./sqrt(obs(1).weight);
tmp.err(2:3:3*nobs) = obs(1).nerr./sqrt(obs(1).weight);
tmp.err(3:3:3*nobs) = obs(1).herr;
%
d(1).ind = find(tmp.err ~= 0)';
d(1).obs = tmp.obs(d(1).ind)';
d(1).err = tmp.err(d(1).ind)';
d(1).mid = [];
d(1).cnt = 0;
%
% (G(1).C * (( G(1).T * ( G(1).B1 - G(1).B2 ) * Mp)*Mc ) + G(1).P * Mp
% 
% (G(1).C * (( G(1).T * ( G(1).B1 - G(1).B2 ) * Mp)*Mc ) + G(1).P * Mp +
% G(1).I * Mi  % including internal deformation
%
G(1).t  = zeros(3*blk(1).nb,2.*blk(1).nb);
G(1).b  = zeros(2*blk(1).nb,3.*blk(1).nblock);
G(1).s  = zeros(3*blk(1).nb,3.*blk(1).nb);
tmp.p   = zeros(3*nobs,3.*blk(1).nblock);
tmp.i   = zeros(3*nobs,3.*blk(1).nblock);
tmp.cf  = ones(3*blk(1).nb,1);
tmp.inv = zeros(3*blk(1).nb,1);
% 
mc = 1;
mt = 1;
mr = 1;
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    nf = size(tri(1).bound(nb1,nb2).clon,2);
    if nf ~= 0
      d(1).cnt             = d(1).cnt+1;
      d(1).mid             = zeros(blk(1).nb,1);
      d(1).mid(mr:mr+nf-1) = 1;
      d(1).mid             = [d(1).mid d(1).mid];
      tmp.c(1:3*nobs,mc     :mc+  nf-1) = tri(1).bound(nb1,nb2).gustr;
      tmp.c(1:3*nobs,mc+  nf:mc+2*nf-1) = tri(1).bound(nb1,nb2).gudip;
      tmp.c(1:3*nobs,mc+2*nf:mc+3*nf-1) = tri(1).bound(nb1,nb2).gutns;
      tmp.cf( mc   :mc+3*nf-1) = tri(1).bound(nb1,nb2).cf ;
      tmp.inv(mc   :mc+3*nf-1) = tri(1).bound(nb1,nb2).inv;
      G(1).t(mc   :mc+  nf-1,mt   :mt+  nf-1) = diag(tri(1).bound(nb1,nb2).st(:,1));
      G(1).t(mc+nf:mc+2*nf-1,mt   :mt+  nf-1) = diag(tri(1).bound(nb1,nb2).dp(:,1));
      G(1).t(mc   :mc+  nf-1,mt+nf:mt+2*nf-1) = diag(tri(1).bound(nb1,nb2).st(:,2));
      G(1).t(mc+nf:mc+2*nf-1,mt+nf:mt+2*nf-1) = diag(tri(1).bound(nb1,nb2).dp(:,2));
      G(1).b(mt   :mt+  nf-1,3*nb1-2) = -1.*(-tri(1).bound(nb1,nb2).oxyz(:,7).*tri(1).bound(nb1,nb2).oxyz(:,3));
      G(1).b(mt   :mt+  nf-1,3*nb1-1) = -1.*(-tri(1).bound(nb1,nb2).oxyz(:,5).*tri(1).bound(nb1,nb2).oxyz(:,3));
      G(1).b(mt   :mt+  nf-1,3*nb1  ) = -1.*( tri(1).bound(nb1,nb2).oxyz(:,5).*tri(1).bound(nb1,nb2).oxyz(:,2)...
                                      +tri(1).bound(nb1,nb2).oxyz(:,7).*tri(1).bound(nb1,nb2).oxyz(:,1));
      G(1).b(mt+nf:mt+2*nf-1,3*nb1-2) = -1.*( tri(1).bound(nb1,nb2).oxyz(:,4).*tri(1).bound(nb1,nb2).oxyz(:,5).*tri(1).bound(nb1,nb2).oxyz(:,3)...
                                      +tri(1).bound(nb1,nb2).oxyz(:,6).*tri(1).bound(nb1,nb2).oxyz(:,2));
      G(1).b(mt+nf:mt+2*nf-1,3*nb1-1) = -1.*(-tri(1).bound(nb1,nb2).oxyz(:,4).*tri(1).bound(nb1,nb2).oxyz(:,7).*tri(1).bound(nb1,nb2).oxyz(:,3)...
                                      -tri(1).bound(nb1,nb2).oxyz(:,6).*tri(1).bound(nb1,nb2).oxyz(:,1));
      G(1).b(mt+nf:mt+2*nf-1,3*nb1  ) = -1.*( tri(1).bound(nb1,nb2).oxyz(:,4).*tri(1).bound(nb1,nb2).oxyz(:,7).*tri(1).bound(nb1,nb2).oxyz(:,2)...
                                      -tri(1).bound(nb1,nb2).oxyz(:,4).*tri(1).bound(nb1,nb2).oxyz(:,5).*tri(1).bound(nb1,nb2).oxyz(:,1));
      G(1).b(mt   :mt+  nf-1,3*nb2-2) = -G(1).b(mt   :mt+  nf-1,3*nb1-2);
      G(1).b(mt   :mt+  nf-1,3*nb2-1) = -G(1).b(mt   :mt+  nf-1,3*nb1-1);
      G(1).b(mt   :mt+  nf-1,3*nb2  ) = -G(1).b(mt   :mt+  nf-1,3*nb1  );
      G(1).b(mt+nf:mt+2*nf-1,3*nb2-2) = -G(1).b(mt+nf:mt+2*nf-1,3*nb1-2);
      G(1).b(mt+nf:mt+2*nf-1,3*nb2-1) = -G(1).b(mt+nf:mt+2*nf-1,3*nb1-1);
      G(1).b(mt+nf:mt+2*nf-1,3*nb2  ) = -G(1).b(mt+nf:mt+2*nf-1,3*nb1  );
      G(1).s(mc     :mc+  nf-1,mc     :mc+  nf-1) = tri(1).bound(nb1,nb2).gsstrT(1:3:end,:);
      G(1).s(mc     :mc+  nf-1,mc+  nf:mc+2*nf-1) = tri(1).bound(nb1,nb2).gsdipT(1:3:end,:);
      G(1).s(mc     :mc+  nf-1,mc+2*nf:mc+3*nf-1) = tri(1).bound(nb1,nb2).gstnsT(1:3:end,:);
      G(1).s(mc+  nf:mc+2*nf-1,mc     :mc+  nf-1) = tri(1).bound(nb1,nb2).gsstrT(2:3:end,:);
      G(1).s(mc+  nf:mc+2*nf-1,mc+  nf:mc+2*nf-1) = tri(1).bound(nb1,nb2).gsdipT(2:3:end,:);
      G(1).s(mc+  nf:mc+2*nf-1,mc+2*nf:mc+3*nf-1) = tri(1).bound(nb1,nb2).gstnsT(2:3:end,:);
      G(1).s(mc+2*nf:mc+3*nf-1,mc     :mc+  nf-1) = tri(1).bound(nb1,nb2).gsstrT(3:3:end,:);
      G(1).s(mc+2*nf:mc+3*nf-1,mc+  nf:mc+2*nf-1) = tri(1).bound(nb1,nb2).gsdipT(3:3:end,:);
      G(1).s(mc+2*nf:mc+3*nf-1,mc+2*nf:mc+3*nf-1) = tri(1).bound(nb1,nb2).gstnsT(3:3:end,:);
%       
      mc = mc + 3*nf;
      mt = mt + 2*nf;
      mr = mr +   nf;
    end
  end
%   
  ind  = obs(1).ablk == nb1;  
  nind = [zeros(size(ind)),ind,zeros(size(ind))]; nind=logical(reshape(nind',3*nobs,1));
  eind = [ind,zeros(size(ind)),zeros(size(ind))]; eind=logical(reshape(eind',3*nobs,1));
  tmp.p(eind,3*nb1-2) = -obs(1).axyz(ind,7).*obs(1).axyz(ind,3);
  tmp.p(eind,3*nb1-1) = -obs(1).axyz(ind,5).*obs(1).axyz(ind,3);
  tmp.p(eind,3*nb1  ) =  obs(1).axyz(ind,5).*obs(1).axyz(ind,2)                    +obs(1).axyz(ind,7).*obs(1).axyz(ind,1);
  tmp.p(nind,3*nb1-2) =  obs(1).axyz(ind,4).*obs(1).axyz(ind,5).*obs(1).axyz(ind,3)+obs(1).axyz(ind,6).*obs(1).axyz(ind,2);
  tmp.p(nind,3*nb1-1) = -obs(1).axyz(ind,4).*obs(1).axyz(ind,7).*obs(1).axyz(ind,3)-obs(1).axyz(ind,6).*obs(1).axyz(ind,1);
  tmp.p(nind,3*nb1  ) =  obs(1).axyz(ind,4).*obs(1).axyz(ind,7).*obs(1).axyz(ind,2)-obs(1).axyz(ind,4).*obs(1).axyz(ind,5).*obs(1).axyz(ind,1);
  tmp.i(eind,3*nb1-2) =  (blk(nb1).xinter).*10^6;
  tmp.i(eind,3*nb1-1) =  (blk(nb1).yinter).*10^6;
  tmp.i(eind,3*nb1  ) =  0;
  tmp.i(nind,3*nb1-2) =  0;
  tmp.i(nind,3*nb1-1) =  (blk(nb1).xinter).*10^6;
  tmp.i(nind,3*nb1  ) =  (blk(nb1).yinter).*10^6;
end
% 
G(1).c     = tmp.c(d(1).ind,:);
G(1).p     = tmp.p(d(1).ind,:);
G(1).i     = tmp.i(d(1).ind,:);
G(1).tb    = sparse(G(1).t*G(1).b);
d(1).mid   = logical(repmat(d(1).mid,3,1));
d(1).cfinv = tmp.cf.*tmp.inv;
end

%% Define initial locking patches
function [d] = InitialLockingPatch(blk,tri,d)
% Coded by Hiroshi Kimura in 2019/2/14
% d(1).id_lock : Locking, give back-slip  (coupling = 1)
% d(1).id_crep : Creeping, slip passively (coupling = 0)
% 
d(1).id_lock = zeros(3*tri(1).nb,1);
d(1).id_crep = zeros(3*tri(1).nb,1);
mc = 1;
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    nf = size(tri(1).bound(nb1,nb2).clon,2);
    if nf ~= 0
      for np = 1:size(blk(1).bound(nb1,nb2).patch,2)
        slipid = inpolygon(tri(1).bound(nb1,nb2).clon,...
                           tri(1).bound(nb1,nb2).clat,...
                           blk(1).bound(nb1,nb2).patch(np).lon,...
                           blk(1).bound(nb1,nb2).patch(np).lat);
        d(1).id_lock(mc:mc+nf-1) = d(1).id_lock(mc:mc+nf-1) | slipid';
%         d(1).id_lock(mc:mc+3*nf-1) = d(1).id_lock(mc:mc+3*nf-1) | [repmat( slipid',2,1); false(nf,1)];
%         d(1).id_crep(mc:mc+3*nf-1) = d(1).id_crep(mc:mc+3*nf-1) | [repmat(~slipid',2,1); false(nf,1)];
      end
      d(1).id_lock(mc:mc+3*nf-1) = [repmat( d(1).id_lock(mc:mc+nf-1),2,1); false(nf,1)];
      d(1).id_crep(mc:mc+3*nf-1) = [repmat(~d(1).id_lock(mc:mc+nf-1),2,1); false(nf,1)];
      mc = mc + 3*nf;
    end
  end
end
d(1).id_lock = logical(d(1).id_lock);
d(1).id_crep = logical(d(1).id_crep);

end

%% Calculate passive slip by locked patch
function [cal] = CalcSlip(blk,tri,prm,obs,eul,d,G)
% Test version coded by Hiroshi Kimura in 2019/2/1
precision  = 'double'     ;
rwd        = prm.rwd      ;
nb         = blk(1).nblock;
passivelim = 1            ;  % [mm], TODO : How to determine?

% initial parameters
mc.int = 1e-2 ;
mp.int = 1e-10;
mi.int = 1e-10;
mc.n   = blk(1).nb;
mp.n   = 3.*blk(1).nblock;
mi.n   = 3.*blk(1).nblock;
% substitute euler pole vectors
mp.old         = double(blk(1).pole);
mp.old(eul.id) = 0                  ;
mp.old         = mp.old+eul.fixw    ;
% substitute internal strain tensors
mi.old = 1e-10.*(-0.5+rand(mi.n,1,precision));
mi.old = mi.old.*blk(1).idinter              ;

% substitude index of locking patches
id_lock = logical(d(1).id_lock);
id_crep = logical(d(1).id_crep);

% Calculate back-slip on locked patches.
cal.slip             = (G(1).tb*mp.old).*d(1).cfinv.*id_lock;

count = 1;
while 1
  % Calculate strain out of locked patches.
  cal.strain           = (G(1).s*cal.slip).*id_crep;
  % Inverse velocity out of locked patches.
  cal.slip(id_crep) = -(G(1).s(id_crep,id_crep) \ cal.strain(id_crep));

  % Sampling initial stain and slip
  if count == 1
    cal.intstrain = cal.strain;
    cal.intslip   = cal.slip  ;
  end
  
  % Find back-slipped region
  id_lock_next = zeros(3*tri(1).nb,1);
  id_crep_next = zeros(3*tri(1).nb,1);
  mc = 1;
  for nb1 = 1:blk(1).nblock
    for nb2 = nb1+1:blk(1).nblock
      nf = size(tri(1).bound(nb1,nb2).clon,2);
      if nf ~= 0
        slipid = sqrt(cal.slip(mc:mc+nf-1).^2+cal.slip(mc+nf:mc+2*nf-1).^2) >= passivelim;
        id_lock_next(mc:mc+3*nf-1) = id_lock_next(mc:mc+3*nf-1) | [repmat( slipid,2,1); false(nf,1)];
        id_crep_next(mc:mc+3*nf-1) = id_crep_next(mc:mc+3*nf-1) | [repmat(~slipid,2,1); false(nf,1)];
        mc = mc + 3*nf;
      end
    end
  end
  id_lock_next = logical(id_lock_next);
  id_crep_next = logical(id_crep_next);
  %
  if sum(id_lock_next==id_lock)==3*tri(1).nb  % Passive slip is almost 0
    break
  else
    count   = count + 1   ;
    id_lock = id_lock_next;
    id_crep = id_crep_next;    
  end
  %
end
cal.it_count = count;
% Rigid motion
cal.rig = G.p*mp.old;
% Elastic motion due to slip deficit
cal.ela = G.c*cal.slip;
% Internal deformation
cal.ine = G.i*mi.old;
% Total velocity
cal.smp = cal.rig+cal.ela+cal.ine;
% Make figure
MakeFigs(blk,tri,cal,obs);

% 
%{  
%  TO DO
%  MCMC inversion
%  calculate velocities
cal.rig = G.p*mp.smp;
cal.ela = G.c*((G.tb*mp.smp).*d(1).cfinv.*mc.smpmat);
cal.ine = G.i*mi.smp;
cal.smp = cal.rig+cal.ela+cal.ine;

mc.std=mc.int.*ones(mc.n,1,precision);
mp.std=mp.int.*ones(mp.n,1,precision);
mi.std=mi.int.*ones(mi.n,1,precision);
mp.old= double(blk(1).pole);
mi.old= 1e-10.*(-0.5+rand(mi.n,1,precision));
cha.mc= zeros(mc.n,prm.kep,precision);
cha.mp= zeros(mp.n,prm.kep,precision);
cha.mi= zeros(mi.n,prm.kep,precision);
% set fix poles if pol.fixflag=1
mp.old(pol.id)=0; mp.old=mp.old+pol.fixw;
mp.std(pol.id)=0;
%
mi.old=mi.old.*blk(1).idinter;
mi.std=mi.std.*blk(1).idinter;
% 
res.old=inf(1,1,precision);
rwdscale=1000*rwd/(prm.cha);
mcscale=rwdscale*0.13;
mpscale=rwdscale*(1.3e-9).*ones(mp.n,1,precision).*~pol.id;
miscale=rwdscale*1e-10;

mc.smpmat=repmat(mc.smp,3,d.cnt);
mc.smpmat=mc.smpmat(d.mid);
%}

end

%% Estimate mechanical coupled area by MCMC (Metropolis-Hasting)
function MechCoupling_MCMC_MH
% Initial value

% Generate next sample

% Accept or reject

end

%% Make figures
function MakeFigs(blk,tri,cal,obs)
% 
figure(100); clf(100)
mc = 1;
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    nf = size(tri(1).bound(nb1,nb2).clon,2);
    if nf ~= 0
      % Back-slip rate
      patch(blk(1).bound(nb1,nb2).blon',...
            blk(1).bound(nb1,nb2).blat',...
      sqrt(cal.slip(mc:mc+nf-1).^2+cal.slip(mc+nf:mc+2*nf-1).^2))
      % Initial patches
      for np = 1:size(blk(1).bound(nb1,nb2).patch,2)
        hold on
        plot(blk(1).bound(nb1,nb2).patch(np).lon,...
             blk(1).bound(nb1,nb2).patch(np).lat,...
             'LineWidth',3,'Color','b')
      end
      mc = mc + 3*nf;
    end
  end
end
colormap('hot')
colorbar
% Velocity at surface due to back-slip
hold on
quiver(obs(1).alon,obs(1).alat,obs(1).evec,obs(1).nvec,'green')
hold on
quiver(obs(1).alon,obs(1).alat,cal.smp(1:3:end)',cal.smp(2:3:end)','blue')
% 
end

%% PLTXY
%====================================================
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
%% XYTPL
%====================================================
function [LAT,LON]=XYTPL(X,Y,ALAT0,ALON0)
%-------------------------------------------------------------------------
%  PLTXY TRANSFORMS (X,Y) TO (ALAT,ALONG)
%  TRANSFORMATION  BETWEEN (X,Y) AND (ALAT,ALONG).
%-------------------------------------------------------------------------
A=6.378160e3;
E2=6.6944541e-3;
E12=6.7395719e-3;
D=5.72958e1;
RD=1.0/D;
RLATO = ALAT0.*RD;
SLATO = sin(RLATO);
R     = A.*(1-E2)./sqrt((1-E2.*SLATO.^2).^3);
AN    = A./sqrt(1.0-E2.*SLATO.^2);
V2    = 1 + E12.*cos(RLATO).^2;
C1    = D./R;
C2    = D./AN;
PH1   = ALAT0+C1.*Y;
RPH1  = PH1.*RD;
TPHI1 = tan(RPH1);
CPHI1 = cos(RPH1);
LAT   = PH1-(C2.*X).^2.*V2.*TPHI1./(2.*D);
LON   = ALON0+C2.*X./CPHI1-(C2.*X).^3.*(1.0+2.*TPHI1.^2)./(6.*D.^2.*CPHI1);
end

%% CONVERT TO XYZ FROM ELL
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

%% CONVERT TO XYZ FROM ELL AT SURFACE
function [OOxyz]=conv2ell(Olat,Olon,Ohig)
Olat=Olat(:);
Olon=Olon(:);
Ohig=Ohig(:);
deg2rad=pi/180;
[Oxyz(:,1),Oxyz(:,2),Oxyz(:,3)]=ell2xyz(Olat,Olon,Ohig);
Oxyz = Oxyz*1e3;
OOxyz=[Oxyz sin(Olat*deg2rad) sin(Olon*deg2rad) cos(Olat*deg2rad) cos(Olon*deg2rad)];
end

%% Convert from (x,y,z) to (lon,lat,ang)
function [lat,lon,ang]=xyzp2lla(X,Y,Z)
% XYZP2LLA  Converts Shpear coordinates from cartesian. Vectorized.
% GRS80
% CODE BY T.ITO 2017/03/11     ver0.1
% lat: deg, lon: deg, ang: deg/m.y.
lat=atan2(Z,sqrt(X.*X+Y.*Y)).*180/pi;
lon=atan2(Y,X).*180/pi;
ang=sqrt(X.*X+Y.*Y+Z.*Z).*(1e6.*(180./pi));
end
