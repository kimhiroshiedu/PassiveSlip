%% PassiveSlip.m
function PassiveSlipKim
% Coded by Ryohei Sasajima final 2013/12/23
% Combined by Hiroshi Kimura 2018/11/12
%--- test
prm.input = 'PARAMETER/parameter_test.txt';
prm.optfile='PARAMETER/opt_bound_par.txt';
prm.alon0 = 0;
prm.alat0 = 0;
%--

% [triC,tri3,tri,sll]=make_test_trill;
% 
% ALON0=142.5;
% ALAT0=38.3;
% [trixyzC,trixyz3,sxyz]=trill2trixyz(triC,tri3,sll,ALAT0,ALON0);
% 
% [sitaS,sitaD,normVec]=strike_dip(trixyzC,trixyz3);

% [xyz]=makexyz;
[prm]     = ReadParameters(prm);
[obs]     = ReadObs(prm);
[blk,obs] = ReadBlockBound(prm,obs);
[blk]     = ReadBlockInterface(blk,prm);
[eul,prm] = ReadEulerPoles(blk,prm);
[blk,prm] = ReadDippingBound(blk,prm);
[blk]     = ReadLockedPatch(blk,prm);

[tri]     = MakeGreenFunction(blk,obs);

[Gu] = makeGreenDisp(obs,trixyz3);
[Gs] = makeGreenStrain(trixyzC,trixyz3,sitaS,sitaD,normVec);
sUxyz=Gu.st;
dUxyz=Gu.dp;
sSsn=Gs.stst;
sSdn=Gs.stdp;
dSsn=Gs.dpst;
dSdn=Gs.dpdp;

[Slip]=defineSlipQ(triC,sitaS);

[xySlip]=Slip2xyll(Slip,sitaS);
figure(103); quiver(triC(:,1),triC(:,2),xySlip(:,1),-xySlip(:,2),'b')
figure(101); quiver(trixyzC(:,1),-trixyzC(:,2),xySlip(:,1),-xySlip(:,2),'g')

ll(:,1:2)=triC(:,1:2);
[vnorm,v,NARFA]=SasaEular2Velo(ll);

%{
[]=relavec2vecVxy(NARFA)
[xyBslip]=BackSlipxy(mcmcC,triC,trixyzC,sitaS,accuV,vecVxy)
%}

[A_Ssn,A_Sdn]=I_CalcStrain(Slip,sSsn,sSdn,dSsn,dSdn);
figure(201); quiver(trixyzC(:,1),-trixyzC(:,2),A_Ssn,A_Sdn,5,'b')

[F_Slip]=OutofAsperity01(A_Ssn,A_Sdn,Slip,sSsn,sSdn,dSsn,dSdn);

[xyFSlip]=FSlip2xyll(F_Slip,sitaS);

Fid=fopen('./QxySlip11.dat','w');
n=size(Slip,1);

for i=1:n
    fprintf(Fid, '%11.3f %11.3f %9.5f %9.5f\n',trixyzC(i,1), trixyzC(i,2), xyFSlip(i,1), xyFSlip(i,2));
end
fclose(Fid);

figure(1); quiver(trixyzC(:,1),-trixyzC(:,2),xyFSlip(:,1),-xyFSlip(:,2),2,'g')
hold on
triplot(tri,1000.*sxyz(:,1),1000.*sxyz(:,2));

[F_Ssn,F_Sdn]=I_CalcStrainF(F_Slip,sSsn,sSdn,dSsn,dSdn);
figure(401); quiver(trixyzC(:,1),-trixyzC(:,2),F_Ssn,F_Sdn,5,'b')

absO_Ssn=abs(O_Ssn);
absO_Sdn=abs(O_Sdn);
sumAOSsn=sum(absO_Ssn,1);
sumAOSdn=sum(absO_Sdn,1);

[FO_Ssn,FO_Sdn]=OutofAsperity11(F_Ssn,F_Sdn,Slip);
absFO_Ssn=abs(FO_Ssn);
absFO_Sdn=abs(FO_Sdn);
sumF0Ssn=sum(absFO_Ssn,1);
sumF0Sdn=sum(absFO_Sdn,1);

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
filerigb           = fscanf(fid,'%s \n',[1,1]); [~] = fgetl(fid);
fileinternal       = fscanf(fid,'%s \n',[1,1]); [~] = fgetl(fid);
dirresult          = fscanf(fid,'%s \n',[1,1]); [~] = fgetl(fid);
prm.home_d = pwd;
prm.fileobs            = fullfile(prm.home_d,fileobs);
prm.dirblock           = fullfile(prm.home_d,dirblock);
prm.dirblock_interface = fullfile(prm.home_d,dirblock_interface);
prm.dirblock_patch     = fullfile(prm.home_d,dirblock_patch);
prm.filepole           = fullfile(prm.home_d,filepole);
prm.filerigb           = fullfile(prm.home_d,filerigb);
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
fprintf('File Rigid boundary       : %s \n',prm.filerigb) 
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
function obs = ReadObs(prm)
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
    pre_tri_f = fullfile(dirblk,['trib_',num2str(nb1),'_',num2str(nb2),'.txt']); 
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
      bo_tri_f = fullfile(dirblk,['tribo_',num2str(nb1),'_',num2str(nb2),'.txt']); 
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
      sub_f = fullfile(dirblk,['b_',num2str(nb1),'_',num2str(nb2),'.txt']);
      fid   = fopen(sub_f,'r');
      if fid >= 0
        fprintf('block interface: %2i  %2i \n',nb1,nb2)
        fprintf('read interface boudary shape file : %s \n',sub_f)
        dep_blk = textscan(fid,'%f%f%f'); fclose(fid);
        dep_blk = cell2mat(dep_blk);
        f    = scatteredinterpolant(dep_blk(:,1),dep_blk(:,2),dep_blk(:,3));
        bo_f = fullfile(dirblk,['bo_',num2str(nb1),'_',num2str(nb2),'.txt']);
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
        out_tri_f = fullfile(prm.dirblock,['trib_',num2str(nb1),'_',num2str(nb2),'.out']);
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
% test version coded by H. Kimura 2019/01/29
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    blk(1).bounc(nb1,nb2).patchid = 0;
    patchfile = fullfile(prm.dirblock_patch,['patchb_',num2str(nb1),'_',num2str(nb2),'.txt']);
    fid       = fopen(patchfile,'r');
    if fid >= 0
      np    = 0;
      tline = fgetl(fid);
      while 1
        if tline < 0
          break
        elseif tline == ">"
          np = np+1;
          blk(1).bound(nb1,nb2).patch(np).lon = [];
          blk(1).bound(nb1,nb2).patch(np).lat = [];
          while 1
            tline = fgetl(fid);
            if tline < 0
              break
            end
            tmp = strtrim(strsplit(tline));
            if ~or(strcmpi(tmp(1),'>'), strcmpi(tmp(1),''))
              tmp = str2double(char(tmp));
              blk(1).bound(nb1,nb2).patch(np).lon = [blk(1).bound(nb1,nb2).patch(np).lon, tmp(1)];
              blk(1).bound(nb1,nb2).patch(np).lat = [blk(1).bound(nb1,nb2).patch(np).lat, tmp(1)];
            else
              break
            end
          end
        elseif tline == ""
          break
        end
      end
    end
  end
end

end
    
%%
function tmp      
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    blk(1).bound(nb1,nb2).type = 1;
    pre_tri_f = fullfile(dirblk,['trib_',num2str(nb1),'_',num2str(nb2),'.txt']); 
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
      bo_tri_f = fullfile(dirblk,['tribo_',num2str(nb1),'_',num2str(nb2),'.txt']); 
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
      sub_f = fullfile(dirblk,['b_',num2str(nb1),'_',num2str(nb2),'.txt']);
      fid   = fopen(sub_f,'r');
      if fid >= 0
        fprintf('block interface: %2i  %2i \n',nb1,nb2)
        fprintf('read interface boudary shape file : %s \n',sub_f)
        dep_blk = textscan(fid,'%f%f%f'); fclose(fid);
        dep_blk = cell2mat(dep_blk);
        f    = scatteredinterpolant(dep_blk(:,1),dep_blk(:,2),dep_blk(:,3));
        bo_f = fullfile(dirblk,['bo_',num2str(nb1),'_',num2str(nb2),'.txt']);
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
        out_tri_f = fullfile(prm.dirblock,['trib_',num2str(nb1),'_',num2str(nb2),'.out']);
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
if exist(prm.FileDipb,'file') ~= 2; return; end  % File not exist
fid = fopen(prm.FileDipb,'r');
tmp = fscanf(fid,'%d %d %d\n',[3 Inf]);
blk(1).dipbo = tmp';
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
    [x(2) x(3)]               = swap(x(2), x(3));
    [y(2) y(3)]               = swap(y(2), y(3));
    [z(2) z(3)]               = swap(z(2), z(3));
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
    [rx ry]                  = RotateXyVec(x(iTri+1)-x(iTri), y(iTri+1)-y(iTri), -strike);
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
        [sx1 sy1]                 = RotateXyVec(sx-x(iTri), sy-y(iTri), -strike);
        sx1(abs(sx1)<0.0001)=0.0001;%for eliminate NAN by R. Sasajima and T. Ito ,Nagoya. U. in 2012%
        sy1(abs(sy1)<0.0001)=0.0001;
        [ux1 uy1 uz1]             = adv(sx1, sy1, sz1-z(iTri), z(iTri), beta, pr, lss, lts, lds);
        
        % Second angular dislocation
        [sx2 sy2]                 = RotateXyVec(sx-x(iTri+1), sy-y(iTri+1), -strike);
        sx2(abs(sx2)<0.0001)=0.0001;%for eliminate NAN by R. Sasajima and T. Ito ,Nagoya. U. in 2012%
        sy2(abs(sy2)<0.0001)=0.0001;
        [ux2 uy2 uz2]             = adv(sx2, sy2, sz1-z(iTri+1), z(iTri+1), beta, pr, lss, lts, lds);
        
        % Rotate vectors to correct for strike
        [uxn uyn]                 = RotateXyVec(ux1-ux2, uy1-uy2, strike);
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
    d                            = LinePlaneIntersect(x, y, z, sx(inPolyIdx(iIdx)), sy(inPolyIdx(iIdx)), sz(inPolyIdx(iIdx)));
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
if denominator == 0;
    denominator                  = eps;
end
t                               = numerator/denominator; % parametric curve parameter
d                               = [sx sy sz]-([sx sy 0]-[sx sy sz])*t;
end

function [xp yp] = RotateXyVec(x, y, alpha)
% Rotate a vector by an angle alpha
x                             = x(:);
y                             = y(:);
alpha                         = pi/180*alpha;
xp                            = cos(alpha).*x - sin(alpha).*y;
yp                            = sin(alpha).*x + cos(alpha).*y;
end

function [v1 v2 v3] = adv(y1, y2, y3, a, beta, nu, B1, B2, B3)
% These are the displacements in a uniform elastic half space due to slip
% on an angular dislocation (Comninou and Dunders, 1975).  Some of the
% equations for the B2 and B3 cases have been corrected following Thomas
% 1993.  The equations are coded in way such that they roughly correspond
% to each line in original text.  Exceptions have been made where it made
% more sense because of grouping symbols.
sinbeta           = sin(beta);
cosbeta           = cos(beta);
cotbeta           = cot(beta);
z1                = y1.*cosbeta - y3.*sinbeta;
z3                = y1.*sinbeta + y3.*cosbeta;
R2                = y1.*y1 + y2.*y2 + y3.*y3;
R                 = sqrt(R2);
y3bar             = y3 + 2.*a;
z1bar             = y1.*cosbeta + y3bar.*sinbeta;
z3bar             = -y1.*sinbeta + y3bar.*cosbeta;
R2bar             = y1.*y1 + y2.*y2 + y3bar.*y3bar;
Rbar              = sqrt(R2bar);
F                 = -atan2(y2, y1) + atan2(y2, z1) + atan2(y2.*R.*sinbeta, y1.*z1+(y2.*y2).*cosbeta);
Fbar              = -atan2(y2, y1) + atan2(y2, z1bar) + atan2(y2.*Rbar.*sinbeta, y1.*z1bar+(y2.*y2).*cosbeta);

% Case I: Burgers vector (B1,0,0)
v1InfB1           = 2.*(1-nu).*(F+Fbar) - y1.*y2.*(1./(R.*(R-y3)) + 1./(Rbar.*(Rbar+y3bar))) - ...
    y2.*cosbeta.*((R.*sinbeta-y1)./(R.*(R-z3)) + (Rbar.*sinbeta-y1)./(Rbar.*(Rbar+z3bar)));
v2InfB1           = (1-2.*nu).*(log(R-y3)+log(Rbar+y3bar) - cosbeta.*(log(R-z3)+log(Rbar+z3bar))) - ...
    y2.*y2.*(1./(R.*(R-y3))+1./(Rbar.*(Rbar+y3bar)) - cosbeta.*(1./(R.*(R-z3))+1./(Rbar.*(Rbar+z3bar))));
v3InfB1           = y2 .* (1./R - 1./Rbar - cosbeta.*((R.*cosbeta-y3)./(R.*(R-z3)) - (Rbar.*cosbeta+y3bar)./(Rbar.*(Rbar+z3bar))));
v1InfB1           = v1InfB1 ./ (8.*pi.*(1-nu));
v2InfB1           = v2InfB1 ./ (8.*pi.*(1-nu));
v3InfB1           = v3InfB1 ./ (8.*pi.*(1-nu));

v1CB1             = -2.*(1-nu).*(1-2.*nu).*Fbar.*(cotbeta.*cotbeta) + (1-2.*nu).*y2./(Rbar+y3bar) .* ((1-2.*nu-a./Rbar).*cotbeta - y1./(Rbar+y3bar).*(nu+a./Rbar)) + ...
    (1-2.*nu).*y2.*cosbeta.*cotbeta./(Rbar+z3bar).*(cosbeta+a./Rbar) + a.*y2.*(y3bar-a).*cotbeta./(Rbar.*Rbar.*Rbar) + ...
    y2.*(y3bar-a)./(Rbar.*(Rbar+y3bar)).*(-(1-2.*nu).*cotbeta + y1./(Rbar+y3bar) .* (2.*nu+a./Rbar) + a.*y1./(Rbar.*Rbar)) + ...
    y2.*(y3bar-a)./(Rbar.*(Rbar+z3bar)).*(cosbeta./(Rbar+z3bar).*((Rbar.*cosbeta+y3bar) .* ((1-2.*nu).*cosbeta-a./Rbar).*cotbeta + 2.*(1-nu).*(Rbar.*sinbeta-y1).*cosbeta) - a.*y3bar.*cosbeta.*cotbeta./(Rbar.*Rbar));
v2CB1             = (1-2.*nu).*((2.*(1-nu).*(cotbeta.*cotbeta)-nu).*log(Rbar+y3bar) -(2.*(1-nu).*(cotbeta.*cotbeta)+1-2.*nu).*cosbeta.*log(Rbar+z3bar)) - ...
    (1-2.*nu)./(Rbar+y3bar).*(y1.*cotbeta.*(1-2.*nu-a./Rbar) + nu.*y3bar - a + (y2.*y2)./(Rbar+y3bar).*(nu+a./Rbar)) - ...
    (1-2.*nu).*z1bar.*cotbeta./(Rbar+z3bar).*(cosbeta+a./Rbar) - a.*y1.*(y3bar-a).*cotbeta./(Rbar.*Rbar.*Rbar) + ...
    (y3bar-a)./(Rbar+y3bar).*(-2.*nu + 1./Rbar.*((1-2.*nu).*y1.*cotbeta-a) + (y2.*y2)./(Rbar.*(Rbar+y3bar)).*(2.*nu+a./Rbar)+a.*(y2.*y2)./(Rbar.*Rbar.*Rbar)) + ...
    (y3bar-a)./(Rbar+z3bar).*((cosbeta.*cosbeta) - 1./Rbar.*((1-2.*nu).*z1bar.*cotbeta+a.*cosbeta) + a.*y3bar.*z1bar.*cotbeta./(Rbar.*Rbar.*Rbar) - 1./(Rbar.*(Rbar+z3bar)) .* ((y2.*y2).*(cosbeta.*cosbeta) - a.*z1bar.*cotbeta./Rbar.*(Rbar.*cosbeta+y3bar)));

v3CB1             = 2.*(1-nu).*(((1-2.*nu).*Fbar.*cotbeta) + (y2./(Rbar+y3bar).*(2.*nu+a./Rbar)) - (y2.*cosbeta./(Rbar+z3bar).*(cosbeta+a./Rbar))) + ...
    y2.*(y3bar-a)./Rbar.*(2.*nu./(Rbar+y3bar)+a./(Rbar.*Rbar)) + ...
    y2.*(y3bar-a).*cosbeta./(Rbar.*(Rbar+z3bar)).*(1-2.*nu-(Rbar.*cosbeta+y3bar)./(Rbar+z3bar).*(cosbeta + a./Rbar) - a.*y3bar./(Rbar.*Rbar));

v1CB1             = v1CB1 ./ (4.*pi.*(1-nu));
v2CB1             = v2CB1 ./ (4.*pi.*(1-nu));
v3CB1             = v3CB1 ./ (4.*pi.*(1-nu));

v1B1              = v1InfB1 + v1CB1;
v2B1              = v2InfB1 + v2CB1;
v3B1              = v3InfB1 + v3CB1;


% Case II: Burgers vector (0,B2,0)
v1InfB2           = -(1-2.*nu).*(log(R-y3) + log(Rbar+y3bar)-cosbeta.*(log(R-z3)+log(Rbar+z3bar))) + ...
    y1.*y1.*(1./(R.*(R-y3))+1./(Rbar.*(Rbar+y3bar))) + z1.*(R.*sinbeta-y1)./(R.*(R-z3)) + z1bar.*(Rbar.*sinbeta-y1)./(Rbar.*(Rbar+z3bar));
v2InfB2           = 2.*(1-nu).*(F+Fbar) + y1.*y2.*(1./(R.*(R-y3))+1./(Rbar.*(Rbar+y3bar))) - y2.*(z1./(R.*(R-z3))+z1bar./(Rbar.*(Rbar+z3bar)));
v3InfB2           = -(1-2.*nu).*sinbeta.*(log(R-z3)-log(Rbar+z3bar)) - y1.*(1./R-1./Rbar) + z1.*(R.*cosbeta-y3)./(R.*(R-z3)) - z1bar.*(Rbar.*cosbeta+y3bar)./(Rbar.*(Rbar+z3bar));
v1InfB2           = v1InfB2 ./ (8.*pi.*(1-nu));
v2InfB2           = v2InfB2 ./ (8.*pi.*(1-nu));
v3InfB2           = v3InfB2 ./ (8.*pi.*(1-nu));

v1CB2             = (1-2.*nu).*((2.*(1-nu).*(cotbeta.*cotbeta)+nu).*log(Rbar+y3bar) - (2.*(1-nu).*(cotbeta.*cotbeta)+1).*cosbeta.*log(Rbar+z3bar)) + ...
    (1-2.*nu)./(Rbar+y3bar).* (-(1-2.*nu).*y1.*cotbeta+nu.*y3bar-a+a.*y1.*cotbeta./Rbar + (y1.*y1)./(Rbar+y3bar).*(nu+a./Rbar)) - ...
    (1-2.*nu).*cotbeta./(Rbar+z3bar).*(z1bar.*cosbeta - a.*(Rbar.*sinbeta-y1)./(Rbar.*cosbeta)) - a.*y1.*(y3bar-a).*cotbeta./(Rbar.*Rbar.*Rbar) + ...
    (y3bar-a)./(Rbar+y3bar).*(2.*nu + 1./Rbar.*((1-2.*nu).*y1.*cotbeta+a) - (y1.*y1)./(Rbar.*(Rbar+y3bar)).*(2.*nu+a./Rbar) - a.*(y1.*y1)./(Rbar.*Rbar.*Rbar)) + ...
    (y3bar-a).*cotbeta./(Rbar+z3bar).*(-cosbeta.*sinbeta+a.*y1.*y3bar./(Rbar.*Rbar.*Rbar.*cosbeta) + (Rbar.*sinbeta-y1)./Rbar.*(2.*(1-nu).*cosbeta - (Rbar.*cosbeta+y3bar)./(Rbar+z3bar).*(1+a./(Rbar.*cosbeta))));
v2CB2             = 2.*(1-nu).*(1-2.*nu).*Fbar.*cotbeta.*cotbeta + (1-2.*nu).*y2./(Rbar+y3bar).*(-(1-2.*nu-a./Rbar).*cotbeta + y1./(Rbar+y3bar).*(nu+a./Rbar)) - ...
    (1-2.*nu).*y2.*cotbeta./(Rbar+z3bar).*(1+a./(Rbar.*cosbeta)) - a.*y2.*(y3bar-a).*cotbeta./(Rbar.*Rbar.*Rbar) + ...
    y2.*(y3bar-a)./(Rbar.*(Rbar+y3bar)).*((1-2.*nu).*cotbeta - 2.*nu.*y1./(Rbar+y3bar) - a.*y1./Rbar.*(1./Rbar+1./(Rbar+y3bar))) + ...
    y2.*(y3bar-a).*cotbeta./(Rbar.*(Rbar+z3bar)).*(-2.*(1-nu).*cosbeta + (Rbar.*cosbeta+y3bar)./(Rbar+z3bar).*(1+a./(Rbar.*cosbeta)) + a.*y3bar./((Rbar.*Rbar).*cosbeta));
v3CB2             = -2.*(1-nu).*(1-2.*nu).*cotbeta .* (log(Rbar+y3bar)-cosbeta.*log(Rbar+z3bar)) - ...
    2.*(1-nu).*y1./(Rbar+y3bar).*(2.*nu+a./Rbar) + 2.*(1-nu).*z1bar./(Rbar+z3bar).*(cosbeta+a./Rbar) + ...
    (y3bar-a)./Rbar.*((1-2.*nu).*cotbeta-2.*nu.*y1./(Rbar+y3bar)-a.*y1./(Rbar.*Rbar)) - ...
    (y3bar-a)./(Rbar+z3bar).*(cosbeta.*sinbeta + (Rbar.*cosbeta+y3bar).*cotbeta./Rbar.*(2.*(1-nu).*cosbeta - (Rbar.*cosbeta+y3bar)./(Rbar+z3bar)) + a./Rbar.*(sinbeta - y3bar.*z1bar./(Rbar.*Rbar) - z1bar.*(Rbar.*cosbeta+y3bar)./(Rbar.*(Rbar+z3bar))));
v1CB2             = v1CB2 ./ (4.*pi.*(1-nu));
v2CB2             = v2CB2 ./ (4.*pi.*(1-nu));
v3CB2             = v3CB2 ./ (4.*pi.*(1-nu));

v1B2              = v1InfB2 + v1CB2;
v2B2              = v2InfB2 + v2CB2;
v3B2              = v3InfB2 + v3CB2;


% Case III: Burgers vector (0,0,B3)
v1InfB3           = y2.*sinbeta.*((R.*sinbeta-y1)./(R.*(R-z3))+(Rbar.*sinbeta-y1)./(Rbar.*(Rbar+z3bar)));
v2InfB3           = (1-2.*nu).*sinbeta.*(log(R-z3)+log(Rbar+z3bar)) - (y2.*y2).*sinbeta.*(1./(R.*(R-z3))+1./(Rbar.*(Rbar+z3bar)));
v3InfB3           = 2.*(1-nu).*(F-Fbar) + y2.*sinbeta.*((R.*cosbeta-y3)./(R.*(R-z3))-(Rbar.*cosbeta+y3bar)./(Rbar.*(Rbar+z3bar)));
v1InfB3           = v1InfB3 ./ (8.*pi.*(1-nu));
v2InfB3           = v2InfB3 ./ (8.*pi.*(1-nu));
v3InfB3           = v3InfB3 ./ (8.*pi.*(1-nu));

v1CB3             = (1-2.*nu).*(y2./(Rbar+y3bar).*(1+a./Rbar) - y2.*cosbeta./(Rbar+z3bar).*(cosbeta+a./Rbar)) - ...
    y2.*(y3bar-a)./Rbar.*(a./(Rbar.*Rbar) + 1./(Rbar+y3bar)) + ...
    y2.*(y3bar-a).*cosbeta./(Rbar.*(Rbar+z3bar)).*((Rbar.*cosbeta+y3bar)./(Rbar+z3bar).*(cosbeta+a./Rbar) + a.*y3bar./(Rbar.*Rbar));
v2CB3             = (1-2.*nu).*(-sinbeta.*log(Rbar+z3bar) - y1./(Rbar+y3bar).*(1+a./Rbar) + z1bar./(Rbar+z3bar).*(cosbeta+a./Rbar)) + ...
    y1.*(y3bar-a)./Rbar.*(a./(Rbar.*Rbar) + 1./(Rbar+y3bar)) - ...
    (y3bar-a)./(Rbar+z3bar).*(sinbeta.*(cosbeta-a./Rbar) + z1bar./Rbar.*(1+a.*y3bar./(Rbar.*Rbar)) - ...
    1./(Rbar.*(Rbar+z3bar)).*((y2.*y2).*cosbeta.*sinbeta - a.*z1bar./Rbar.*(Rbar.*cosbeta+y3bar)));
v3CB3             = 2.*(1-nu).*Fbar + 2.*(1-nu).*(y2.*sinbeta./(Rbar+z3bar).*(cosbeta + a./Rbar)) + ...
    y2.*(y3bar-a).*sinbeta./(Rbar.*(Rbar+z3bar)).*(1 + (Rbar.*cosbeta+y3bar)./(Rbar+z3bar).*(cosbeta+a./Rbar) + a.*y3bar./(Rbar.*Rbar));
v1CB3             = v1CB3 ./ (4.*pi.*(1-nu));
v2CB3             = v2CB3 ./ (4.*pi.*(1-nu));
v3CB3             = v3CB3 ./ (4.*pi.*(1-nu));

v1B3              = v1InfB3 + v1CB3;
v2B3              = v2InfB3 + v2CB3;
v3B3              = v3InfB3 + v3CB3;


% Sum the for each slip component
v1                = B1.*v1B1 + B2.*v1B2 + B3.*v1B3;
v2                = B1.*v2B1 + B2.*v2B2 + B3.*v2B3;
v3                = B1.*v3B1 + B2.*v3B2 + B3.*v3B3;
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
    [x(2) x(3)]               = swap(x(2), x(3));
    [y(2) y(3)]               = swap(y(2), y(3));
    [z(2) z(3)]               = swap(z(2), z(3));
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
    [rx ry]                  = RotateXyVec(x(iTri+1)-x(iTri), y(iTri+1)-y(iTri), -strike);
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
        [sx1 sy1]                 = RotateXyVec(sx-x(iTri), sy-y(iTri), -strike);
        sx1(abs(sx1)<0.0001)=0.0001;%for eliminate NAN by R. Sasajima and T. Ito ,Nagoya. U. in 2012%
        sy1(abs(sy1)<0.0001)=0.0001;
        [a11 a22 a33 a12 a13 a23] = advs(sx1, sy1, sz-z(iTri), z(iTri), beta, pr, lss, lts, lds);
        
        % Second angular dislocation
        [sx2 sy2]                 = RotateXyVec(sx-x(iTri+1), sy-y(iTri+1), -strike);
        sx2(abs(sx2)<0.0001)=0.0001;
        sy2(abs(sy2)<0.0001)=0.0001;
        [b11 b22 b33 b12 b13 b23] = advs(sx2, sy2, sz-z(iTri+1), z(iTri+1), beta, pr, lss, lts, lds);
        
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

function [a b] = swap(a, b)
% Swap two values
temp                            = a;
a                               = b;
b                               = temp;
end

function [e11 e22 e33 e12 e13 e23] = advs(y1, y2, y3, a, b, nu, B1, B2, B3)
% These are the strains in a uniform elastic half space due to slip
% on an angular dislocation.  They were calculated by symbolically
% differentiating the expressions for the displacements (Comninou and
% Dunders, 1975, with typos noted by Thomas 1993) then combining the
% elements of the displacement gradient tensor to form the strain tensor.

e11 = B1.*(1./8.*((2-2.*nu).*(2.*y2./y1.^2./(1+y2.^2./y1.^2)-y2./(y1.*cos(b)-y3.*sin(b)).^2.*cos(b)./(1+y2.^2./(y1.*cos(b)-y3.*sin(b)).^2)+(y2./(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).*y1-y2.*(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2.*(2.*y1.*cos(b)-y3.*sin(b)))./(1+y2.^2.*(y1.^2+y2.^2+y3.^2).*sin(b).^2./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2)-y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*cos(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+(y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*y1-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*(2.*y1.*cos(b)+(y3+2.*a).*sin(b)))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2))-y2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))-y1.*y2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).*y1-1./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).^2.*y1-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*y1-1./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*y1)-y2.*cos(b).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b).*y1-1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y1-((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1-sin(b))+(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).*y1-1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y1-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))))./pi./(1-nu)+1./4.*((-2+2.*nu).*(1-2.*nu).*(y2./y1.^2./(1+y2.^2./y1.^2)-y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*cos(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+(y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*y1-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*(2.*y1.*cos(b)+(y3+2.*a).*sin(b)))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2)).*cot(b).^2-(1-2.*nu).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((1-2.*nu-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)-y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1+(1-2.*nu).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1.*cot(b)-1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2))-(1-2.*nu).*y2.*cos(b).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))-(1-2.*nu).*y2.*cos(b).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1-3.*a.*y2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*y1-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((-1+2.*nu).*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2)).*y1-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((-1+2.*nu).*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2)).*y1+y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)+a./(y1.^2+y2.^2+(y3+2.*a).^2)-2.*a.*y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^2)-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+(2-2.*nu).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1).*cos(b))-a.*(y3+2.*a).*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2)).*y1-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+(2-2.*nu).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1).*cos(b))-a.*(y3+2.*a).*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(-cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+(2-2.*nu).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1).*cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*y1.*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1.*cot(b)+(2-2.*nu).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).*y1-1).*cos(b))+2.*a.*(y3+2.*a).*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y1))./pi./(1-nu))+B2.*(1./8.*((-1+2.*nu).*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1./((y1.^2+y2.^2+y3.^2).^(1./2)-y3)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-cos(b).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1-sin(b))./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))))+2.*y1.*(1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))+y1.^2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).*y1-1./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).^2.*y1-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*y1-1./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*y1)+cos(b).*((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+(y1.*cos(b)-y3.*sin(b)).*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b).*y1-1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-(y1.*cos(b)-y3.*sin(b)).*((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y1-(y1.*cos(b)-y3.*sin(b)).*((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1-sin(b))+cos(b).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+(y1.*cos(b)+(y3+2.*a).*sin(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).*y1-1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y1-(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b)))./pi./(1-nu)+1./4.*((1-2.*nu).*(((2-2.*nu).*cot(b).^2+nu)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-((2-2.*nu).*cot(b).^2+1).*cos(b).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((-1+2.*nu).*y1.*cot(b)+nu.*(y3+2.*a)-a+a.*y1.*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1+(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((-1+2.*nu).*cot(b)+a.*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-a.*y1.^2.*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)+2.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-y1.^3./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.^3./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2))+(1-2.*nu).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*((y1.*cos(b)+(y3+2.*a).*sin(b)).*cos(b)-a.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))-(1-2.*nu).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b).^2-a.*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).*y1-1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)+a.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b).*y1)-a.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)+3.*a.*y1.^2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2)-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((1-2.*nu).*y1.*cot(b)+a)-y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-a.*y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((1-2.*nu).*y1.*cot(b)+a).*y1+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1-2.*nu).*cot(b)-2.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+y1.^3./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+y1.^3./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+y1.^3./(y1.^2+y2.^2+(y3+2.*a).^2).^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a-2.*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1+3.*a.*y1.^3./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2))-(y3+a).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(-cos(b).*sin(b)+a.*y1.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)))).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+(y3+a).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b)-3.*a.*y1.^2.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2)./cos(b)+(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).*y1-1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))).*y1+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b).*y1)))./pi./(1-nu))+B3.*(1./8.*y2.*sin(b).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b).*y1-1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y1-((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1-sin(b))+(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).*y1-1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y1-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b)))./pi./(1-nu)+1./4.*((1-2.*nu).*(-y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1+y2.*cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+y2.*cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1)+y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(a./(y1.^2+y2.^2+(y3+2.*a).^2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)).*y1-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-2.*a./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y1-1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1)-y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*y1-y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1-2.*a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y1))./pi./(1-nu));
e22 = B1.*(1./8.*((1-2.*nu).*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y2./((y1.^2+y2.^2+y3.^2).^(1./2)-y3)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-cos(b).*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y2./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))))-2.*y2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-cos(b).*(1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))))-y2.^2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).*y2-1./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).^2.*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*y2-cos(b).*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y2-1./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*y2)))./pi./(1-nu)+1./4.*((1-2.*nu).*(((2-2.*nu).*cot(b).^2-nu)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-((2-2.*nu).*cot(b).^2+1-2.*nu).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))+(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(y1.*cot(b).*(1-2.*nu-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+nu.*(y3+2.*a)-a+y2.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(a.*y1.*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2+2.*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-y2.^3./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y2.^3./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2))+(1-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2+(1-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2+3.*a.*y2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*y1-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(-2.*nu+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((1-2.*nu).*y1.*cot(b)-a)+y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((1-2.*nu).*y1.*cot(b)-a).*y2+2.*y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-y2.^3./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-y2.^3./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-y2.^3./(y1.^2+y2.^2+(y3+2.*a).^2).^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a+2.*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2-3.*a.*y2.^3./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2))-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b).^2-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((1-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)+a.*cos(b))+a.*(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).^2-a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((1-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)+a.*cos(b)).*y2-3.*a.*(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*y2+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).^2-a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)).*y2+1./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(y2.^2.*cos(b).^2-a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)).*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(2.*y2.*cos(b).^2+a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*y2-a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y2)))./pi./(1-nu))+B2.*(1./8.*((2-2.*nu).*(-2./y1./(1+y2.^2./y1.^2)+1./(y1.*cos(b)-y3.*sin(b))./(1+y2.^2./(y1.*cos(b)-y3.*sin(b)).^2)+((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b))+y2.^2./(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b))-2.*y2.^2.*(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2.*cos(b))./(1+y2.^2.*(y1.^2+y2.^2+y3.^2).*sin(b).^2./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2)+1./(y1.*cos(b)+(y3+2.*a).*sin(b))./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))+y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))-2.*y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*cos(b))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2))+y1.*(1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))+y1.*y2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).*y2-1./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).^2.*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*y2)-(y1.*cos(b)-y3.*sin(b))./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-y2.*(-(y1.*cos(b)-y3.*sin(b))./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y2-(y1.*cos(b)-y3.*sin(b))./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*y2-(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y2-(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*y2))./pi./(1-nu)+1./4.*((2-2.*nu).*(1-2.*nu).*(-1./y1./(1+y2.^2./y1.^2)+1./(y1.*cos(b)+(y3+2.*a).*sin(b))./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))+y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))-2.*y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*cos(b))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2)).*cot(b).^2+(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((-1+2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)))-(1-2.*nu).*y2.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((-1+2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+(1-2.*nu).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2.*cot(b)-y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1)-(1-2.*nu).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))+(1-2.*nu).*y2.^2.*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+(1-2.*nu).*y2.^2.*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b)-a.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)+3.*a.*y2.^2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2)+(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((1-2.*nu).*cot(b)-2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)))-y2.^2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((1-2.*nu).*cot(b)-2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)))-y2.^2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((1-2.*nu).*cot(b)-2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)))+y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2+a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)).*y2-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2-1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2))+(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*((-2+2.*nu).*cos(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)./cos(b))-y2.^2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*((-2+2.*nu).*cos(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)./cos(b))-y2.^2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*((-2+2.*nu).*cos(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)./cos(b))+y2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b).*y2-2.*a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^2./cos(b).*y2))./pi./(1-nu))+B3.*(1./8.*((1-2.*nu).*sin(b).*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y2./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-2.*y2.*sin(b).*(1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-y2.^2.*sin(b).*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y2-1./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*y2))./pi./(1-nu)+1./4.*((1-2.*nu).*(-sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1+y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1-(y1.*cos(b)+(y3+2.*a).*sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-(y1.*cos(b)+(y3+2.*a).*sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2)-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(a./(y1.^2+y2.^2+(y3+2.*a).^2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)).*y1+y1.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-2.*a./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y2-1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2)+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(sin(b).*(cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).*sin(b)-a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(sin(b).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2-(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(1+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*y2-2.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*a.*(y3+2.*a).*y2+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).*sin(b)-a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)).*y2+1./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(y2.^2.*cos(b).*sin(b)-a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)).*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(2.*y2.*cos(b).*sin(b)+a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*y2-a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y2)))./pi./(1-nu));
e33 = B1.*(1./8.*y2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2).*y3+1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)-cos(b).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*cos(b).*y3-1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y3-((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-cos(b))-(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+1./2.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(2.*y3+4.*a)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))))./pi./(1-nu)+1./4.*((2-2.*nu).*((1-2.*nu).*(-y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*sin(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+(1./2.*y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*(2.*y3+4.*a)-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*y1)./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2)).*cot(b)-y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)-1./2.*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)+y2.*cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+1./2.*y2.*cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a))+y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*nu./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)+a./(y1.^2+y2.^2+(y3+2.*a).^2))-1./2.*y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*nu./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)+a./(y1.^2+y2.^2+(y3+2.*a).^2)).*(2.*y3+4.*a)+y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-2.*nu./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*(2.*y3+4.*a))+y2.*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1-2.*nu-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-1./2.*y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1-2.*nu-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*(2.*y3+4.*a)-y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1-2.*nu-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(-(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+1./2.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)-a./(y1.^2+y2.^2+(y3+2.*a).^2)+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*(2.*y3+4.*a)))./pi./(1-nu))+B2.*(1./8.*((-1+2.*nu).*sin(b).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-cos(b))./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-y1.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2).*y3+1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a))-sin(b).*((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+(y1.*cos(b)-y3.*sin(b)).*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*cos(b).*y3-1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-(y1.*cos(b)-y3.*sin(b)).*((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y3-(y1.*cos(b)-y3.*sin(b)).*((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-cos(b))-sin(b).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-(y1.*cos(b)+(y3+2.*a).*sin(b)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+1./2.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(2.*y3+4.*a)+(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b)))./pi./(1-nu)+1./4.*((-2+2.*nu).*(1-2.*nu).*cot(b).*((1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-cos(b).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))+(2-2.*nu).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)+1./2.*(2-2.*nu).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)+(2-2.*nu).*sin(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-(2-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))-1./2.*(2-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((1-2.*nu).*cot(b)-2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2))-1./2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((1-2.*nu).*cot(b)-2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2)).*(2.*y3+4.*a)+(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)+a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*(2.*y3+4.*a))-1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b).*sin(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(sin(b)-(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2)-(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))))+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b).*sin(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(sin(b)-(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2)-(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*((1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-1./2.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))).*(2.*y3+4.*a)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b)))-1./2.*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(sin(b)-(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2)-(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))).*(2.*y3+4.*a)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2)-(y3+2.*a).*sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2)+(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*(2.*y3+4.*a)-sin(b).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-(y1.*cos(b)+(y3+2.*a).*sin(b)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+1./2.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(2.*y3+4.*a)+(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b)))))./pi./(1-nu))+B3.*(1./8.*((2-2.*nu).*(y2./(y1.*cos(b)-y3.*sin(b)).^2.*sin(b)./(1+y2.^2./(y1.*cos(b)-y3.*sin(b)).^2)+(y2./(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).*y3+y2.*(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b).^2./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2.*y1)./(1+y2.^2.*(y1.^2+y2.^2+y3.^2).*sin(b).^2./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2)+y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*sin(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)-(1./2.*y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*(2.*y3+4.*a)-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*y1)./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2))+y2.*sin(b).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*cos(b).*y3-1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y3-((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-cos(b))-(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+1./2.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(2.*y3+4.*a)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))))./pi./(1-nu)+1./4.*((2-2.*nu).*(-y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*sin(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+(1./2.*y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*(2.*y3+4.*a)-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*y1)./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2))-(2-2.*nu).*y2.*sin(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))-1./2.*(2-2.*nu).*y2.*sin(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)+y2.*sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-1./2.*y2.*(y3+a).*sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*(2.*y3+4.*a)-y2.*(y3+a).*sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+y2.*(y3+a).*sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*((1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))-1./2.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)+a./(y1.^2+y2.^2+(y3+2.*a).^2)-a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*(2.*y3+4.*a)))./pi./(1-nu));
e12 = 1./2.*B1.*(1./8.*((2-2.*nu).*(-2./y1./(1+y2.^2./y1.^2)+1./(y1.*cos(b)-y3.*sin(b))./(1+y2.^2./(y1.*cos(b)-y3.*sin(b)).^2)+((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b))+y2.^2./(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b))-2.*y2.^2.*(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2.*cos(b))./(1+y2.^2.*(y1.^2+y2.^2+y3.^2).*sin(b).^2./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2)+1./(y1.*cos(b)+(y3+2.*a).*sin(b))./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))+y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))-2.*y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*cos(b))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2))-y1.*(1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))-y1.*y2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).*y2-1./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).^2.*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*y2)-cos(b).*(((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-y2.*cos(b).*(1./(y1.^2+y2.^2+y3.^2).*sin(b).*y2./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y2-((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*y2+1./(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y2-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*y2))./pi./(1-nu)+1./4.*((-2+2.*nu).*(1-2.*nu).*(-1./y1./(1+y2.^2./y1.^2)+1./(y1.*cos(b)+(y3+2.*a).*sin(b))./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))+y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))-2.*y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*cos(b))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2)).*cot(b).^2+(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((1-2.*nu-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)-y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)))-(1-2.*nu).*y2.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((1-2.*nu-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)-y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+(1-2.*nu).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2.*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2+y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1)+(1-2.*nu).*cos(b).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-(1-2.*nu).*y2.^2.*cos(b).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-(1-2.*nu).*y2.^2.*cos(b).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)+a.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)-3.*a.*y2.^2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2)+(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((-1+2.*nu).*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2))-y2.^2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((-1+2.*nu).*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2))-y2.^2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((-1+2.*nu).*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2))+y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1-2.*a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y2)+(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+(2-2.*nu).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1).*cos(b))-a.*(y3+2.*a).*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2))-y2.^2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+(2-2.*nu).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1).*cos(b))-a.*(y3+2.*a).*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2))-y2.^2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+(2-2.*nu).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1).*cos(b))-a.*(y3+2.*a).*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2))+y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(-cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+(2-2.*nu).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1).*cos(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2+cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*y2.*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2.*cot(b)+(2-2.*nu)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).*y2.*cos(b))+2.*a.*(y3+2.*a).*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y2))./pi./(1-nu))+1./2.*B2.*(1./8.*((-1+2.*nu).*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y2./((y1.^2+y2.^2+y3.^2).^(1./2)-y3)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-cos(b).*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y2./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))))+y1.^2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).*y2-1./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).^2.*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*y2)+(y1.*cos(b)-y3.*sin(b))./(y1.^2+y2.^2+y3.^2).*sin(b).*y2./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-(y1.*cos(b)-y3.*sin(b)).*((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y2-(y1.*cos(b)-y3.*sin(b)).*((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*y2+(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y2-(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*y2)./pi./(1-nu)+1./4.*((1-2.*nu).*(((2-2.*nu).*cot(b).^2+nu)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-((2-2.*nu).*cot(b).^2+1).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((-1+2.*nu).*y1.*cot(b)+nu.*(y3+2.*a)-a+a.*y1.*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2+(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-a.*y1.*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2-y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2)+(1-2.*nu).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*((y1.*cos(b)+(y3+2.*a).*sin(b)).*cos(b)-a.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-(1-2.*nu).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(-a./(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).*y2./cos(b)+a.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b).*y2)+3.*a.*y2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*y1-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((1-2.*nu).*y1.*cot(b)+a)-y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-a.*y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((1-2.*nu).*y1.*cot(b)+a).*y2+y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*y2+y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*y2+y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a.*y2+3.*a.*y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*y2)-(y3+a).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(-cos(b).*sin(b)+a.*y1.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2+(y3+a).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(-3.*a.*y1.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2)./cos(b).*y2+1./(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).*y2.*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))).*y2+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b).*y2)))./pi./(1-nu))+1./2.*B3.*(1./8.*sin(b).*(((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))./pi./(1-nu)+1./8.*y2.*sin(b).*(1./(y1.^2+y2.^2+y3.^2).*sin(b).*y2./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y2-((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*y2+1./(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y2-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*y2)./pi./(1-nu)+1./4.*((1-2.*nu).*(1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-y2.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y2.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)-cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+y2.^2.*cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y2.^2.*cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2))-(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(a./(y1.^2+y2.^2+(y3+2.*a).^2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))+y2.^2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(a./(y1.^2+y2.^2+(y3+2.*a).^2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-2.*a./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y2-1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2)+(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-y2.^2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-y2.^2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))+y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2-2.*a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y2))./pi./(1-nu))+1./2.*B1.*(1./8.*((1-2.*nu).*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1./((y1.^2+y2.^2+y3.^2).^(1./2)-y3)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-cos(b).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1-sin(b))./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))))-y2.^2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).*y1-1./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).^2.*y1-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*y1-1./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*y1-cos(b).*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y1-1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1-sin(b))-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y1-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b)))))./pi./(1-nu)+1./4.*((1-2.*nu).*(((2-2.*nu).*cot(b).^2-nu)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-((2-2.*nu).*cot(b).^2+1-2.*nu).*cos(b).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))+(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(y1.*cot(b).*(1-2.*nu-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+nu.*(y3+2.*a)-a+y2.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((1-2.*nu-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+a.*y1.^2.*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)-y2.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-y2.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1)-(1-2.*nu).*cos(b).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+(1-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+(1-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1-a.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)+3.*a.*y1.^2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2)-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(-2.*nu+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((1-2.*nu).*y1.*cot(b)-a)+y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((1-2.*nu).*y1.*cot(b)-a).*y1+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1-2.*nu).*cot(b)-y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*y1-y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*y1-y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a.*y1-3.*a.*y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*y1)-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b).^2-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((1-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)+a.*cos(b))+a.*(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).^2-a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a))).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((1-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)+a.*cos(b)).*y1-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1-2.*nu).*cos(b).*cot(b)+a.*(y3+2.*a).*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)-3.*a.*(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*y1+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).^2-a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)).*y1+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(y2.^2.*cos(b).^2-a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(-a.*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)+a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*y1-a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y1)))./pi./(1-nu))+1./2.*B2.*(1./8.*((2-2.*nu).*(2.*y2./y1.^2./(1+y2.^2./y1.^2)-y2./(y1.*cos(b)-y3.*sin(b)).^2.*cos(b)./(1+y2.^2./(y1.*cos(b)-y3.*sin(b)).^2)+(y2./(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).*y1-y2.*(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2.*(2.*y1.*cos(b)-y3.*sin(b)))./(1+y2.^2.*(y1.^2+y2.^2+y3.^2).*sin(b).^2./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2)-y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*cos(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+(y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*y1-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*(2.*y1.*cos(b)+(y3+2.*a).*sin(b)))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2))+y2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))+y1.*y2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).*y1-1./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).^2.*y1-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*y1-1./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*y1)-y2.*(cos(b)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-(y1.*cos(b)-y3.*sin(b))./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y1-(y1.*cos(b)-y3.*sin(b))./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1-sin(b))+cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y1-(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))))./pi./(1-nu)+1./4.*((2-2.*nu).*(1-2.*nu).*(y2./y1.^2./(1+y2.^2./y1.^2)-y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*cos(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+(y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*y1-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*(2.*y1.*cos(b)+(y3+2.*a).*sin(b)))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2)).*cot(b).^2-(1-2.*nu).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((-1+2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1+(1-2.*nu).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1.*cot(b)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2))+(1-2.*nu).*y2.*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+(1-2.*nu).*y2.*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b).*y1+3.*a.*y2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*y1-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((1-2.*nu).*cot(b)-2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))).*y1-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((1-2.*nu).*cot(b)-2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))).*y1+y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-2.*nu./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)+2.*nu.*y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))+a.*y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1-1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1))-y2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*((-2+2.*nu).*cos(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)./cos(b)).*y1-y2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*((-2+2.*nu).*cos(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)./cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+y2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b).*y1-2.*a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^2./cos(b).*y1))./pi./(1-nu))+1./2.*B3.*(1./8.*((1-2.*nu).*sin(b).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1-sin(b))./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-y2.^2.*sin(b).*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y1-1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1-sin(b))-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y1-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))))./pi./(1-nu)+1./4.*((1-2.*nu).*(-sin(b).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)+cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-(y1.*cos(b)+(y3+2.*a).*sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))-(y1.*cos(b)+(y3+2.*a).*sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1)+(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(a./(y1.^2+y2.^2+(y3+2.*a).^2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))-y1.^2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(a./(y1.^2+y2.^2+(y3+2.*a).^2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))+y1.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-2.*a./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y1-1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1)+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(sin(b).*(cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).*sin(b)-a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a))).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(sin(b).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1+cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(1+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*y1-2.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*a.*(y3+2.*a).*y1+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).*sin(b)-a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)).*y1+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(y2.^2.*cos(b).*sin(b)-a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(-a.*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)+a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*y1-a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y1)))./pi./(1-nu));
e13 = 1./2.*B1.*(1./8.*((2-2.*nu).*(y2./(y1.*cos(b)-y3.*sin(b)).^2.*sin(b)./(1+y2.^2./(y1.*cos(b)-y3.*sin(b)).^2)+(y2./(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).*y3+y2.*(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b).^2./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2.*y1)./(1+y2.^2.*(y1.^2+y2.^2+y3.^2).*sin(b).^2./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2)-y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*sin(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+(1./2.*y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*(2.*y3+4.*a)-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*y1)./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2))-y1.*y2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).*y3-1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-1)-1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*y3+4.*a)-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1))-y2.*cos(b).*(1./(y1.^2+y2.^2+y3.^2).*sin(b).*y3./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y3-((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-cos(b))+1./2./(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).*(2.*y3+4.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-1./2.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(2.*y3+4.*a)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))))./pi./(1-nu)+1./4.*((-2+2.*nu).*(1-2.*nu).*(-y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*sin(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+(1./2.*y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*(2.*y3+4.*a)-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*y1)./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2)).*cot(b).^2-(1-2.*nu).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((1-2.*nu-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)-y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)+(1-2.*nu).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(1./2.*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a).*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)+1./2.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a))-(1-2.*nu).*y2.*cos(b).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))-1./2.*(1-2.*nu).*y2.*cos(b).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2.*cot(b)-3./2.*a.*y2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*(2.*y3+4.*a)+y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((-1+2.*nu).*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2))-1./2.*y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((-1+2.*nu).*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2)).*(2.*y3+4.*a)-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((-1+2.*nu).*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)+y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)-1./2.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*(2.*y3+4.*a))+y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+(2-2.*nu).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1).*cos(b))-a.*(y3+2.*a).*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2))-1./2.*y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+(2-2.*nu).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1).*cos(b))-a.*(y3+2.*a).*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2)).*(2.*y3+4.*a)-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+(2-2.*nu).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1).*cos(b))-a.*(y3+2.*a).*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(-cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+(2-2.*nu).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1).*cos(b)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*((1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1).*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+1./2.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a).*cot(b)+1./2.*(2-2.*nu)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).*(2.*y3+4.*a).*cos(b))-a.*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2)+a.*(y3+2.*a).*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*(2.*y3+4.*a)))./pi./(1-nu))+1./2.*B2.*(1./8.*((-1+2.*nu).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-1)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3)+(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-cos(b).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-cos(b))./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))))+y1.^2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).*y3-1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-1)-1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*y3+4.*a)-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1))-sin(b).*((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+(y1.*cos(b)-y3.*sin(b))./(y1.^2+y2.^2+y3.^2).*sin(b).*y3./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-(y1.*cos(b)-y3.*sin(b)).*((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y3-(y1.*cos(b)-y3.*sin(b)).*((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-cos(b))+sin(b).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+1./2.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).*(2.*y3+4.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-1./2.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(2.*y3+4.*a)-(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b)))./pi./(1-nu)+1./4.*((1-2.*nu).*(((2-2.*nu).*cot(b).^2+nu).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-((2-2.*nu).*cot(b).^2+1).*cos(b).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((-1+2.*nu).*y1.*cot(b)+nu.*(y3+2.*a)-a+a.*y1.*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)+(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu-1./2.*a.*y1.*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)-y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)-1./2.*y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a))+(1-2.*nu).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*((y1.*cos(b)+(y3+2.*a).*sin(b)).*cos(b)-a.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))-(1-2.*nu).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b).*sin(b)-1./2.*a./(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).*(2.*y3+4.*a)./cos(b)+1./2.*a.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b).*(2.*y3+4.*a))-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1.*cot(b)+3./2.*a.*y1.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*(2.*y3+4.*a)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((1-2.*nu).*y1.*cot(b)+a)-y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-a.*y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2))-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((1-2.*nu).*y1.*cot(b)+a)-y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-a.*y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((1-2.*nu).*y1.*cot(b)+a).*(2.*y3+4.*a)+1./2.*y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(2.*y3+4.*a)+y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)+1./2.*y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a.*(2.*y3+4.*a)+3./2.*a.*y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*(2.*y3+4.*a))+cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(-cos(b).*sin(b)+a.*y1.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))))-(y3+a).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(-cos(b).*sin(b)+a.*y1.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)))).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+(y3+a).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b).*y1-3./2.*a.*y1.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2)./cos(b).*(2.*y3+4.*a)+1./2./(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).*(2.*y3+4.*a).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)))-1./2.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))).*(2.*y3+4.*a)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+1./2.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b).*(2.*y3+4.*a))))./pi./(1-nu))+1./2.*B3.*(1./8.*y2.*sin(b).*(1./(y1.^2+y2.^2+y3.^2).*sin(b).*y3./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y3-((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-cos(b))+1./2./(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).*(2.*y3+4.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-1./2.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(2.*y3+4.*a)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b)))./pi./(1-nu)+1./4.*((1-2.*nu).*(-y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)-1./2.*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)+y2.*cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+1./2.*y2.*cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a))-y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(a./(y1.^2+y2.^2+(y3+2.*a).^2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))+1./2.*y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(a./(y1.^2+y2.^2+(y3+2.*a).^2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)).*(2.*y3+4.*a)-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-a./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*(2.*y3+4.*a)-1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1))+y2.*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-1./2.*y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*(2.*y3+4.*a)-y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*((1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))-1./2.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)+a./(y1.^2+y2.^2+(y3+2.*a).^2)-a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*(2.*y3+4.*a)))./pi./(1-nu))+1./2.*B1.*(1./8.*y2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2).*y1+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1-cos(b).*(1./(y1.^2+y2.^2+y3.^2).*cos(b).*y1./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y1-((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1-sin(b))-1./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y1+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))))./pi./(1-nu)+1./4.*((2-2.*nu).*((1-2.*nu).*(y2./y1.^2./(1+y2.^2./y1.^2)-y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*cos(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+(y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*y1-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*(2.*y1.*cos(b)+(y3+2.*a).*sin(b)))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2)).*cot(b)-y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1+y2.*cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+y2.*cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1)-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*nu./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)+a./(y1.^2+y2.^2+(y3+2.*a).^2)).*y1+y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-2.*nu./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-2.*a./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y1)-y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1-2.*nu-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*y1-y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1-2.*nu-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1+2.*a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y1))./pi./(1-nu))+1./2.*B2.*(1./8.*((-1+2.*nu).*sin(b).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1-sin(b))./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-1./(y1.^2+y2.^2+y3.^2).^(1./2)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2).*y1+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1)+cos(b).*((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+(y1.*cos(b)-y3.*sin(b))./(y1.^2+y2.^2+y3.^2).*cos(b).*y1./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-(y1.*cos(b)-y3.*sin(b)).*((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y1-(y1.*cos(b)-y3.*sin(b)).*((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1-sin(b))-cos(b).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y1+(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b)))./pi./(1-nu)+1./4.*((-2+2.*nu).*(1-2.*nu).*cot(b).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-cos(b).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-(2-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+(2-2.*nu).*y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+(2-2.*nu).*y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)+(2-2.*nu).*cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-(2-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))-(2-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1-(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((1-2.*nu).*cot(b)-2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2)).*y1+(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-2.*nu./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)+2.*nu.*y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-a./(y1.^2+y2.^2+(y3+2.*a).^2)+2.*a.*y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^2)+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b).*sin(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(sin(b)-(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2)-(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y1.*cot(b).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))).*y1+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b)))-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(sin(b)-(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2)-(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))).*y1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-(y3+2.*a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2)+2.*(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y1-cos(b).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y1+(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b)))))./pi./(1-nu))+1./2.*B3.*(1./8.*((2-2.*nu).*(-y2./(y1.*cos(b)-y3.*sin(b)).^2.*cos(b)./(1+y2.^2./(y1.*cos(b)-y3.*sin(b)).^2)+(y2./(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).*y1-y2.*(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2.*(2.*y1.*cos(b)-y3.*sin(b)))./(1+y2.^2.*(y1.^2+y2.^2+y3.^2).*sin(b).^2./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2)+y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*cos(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)-(y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*y1-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*(2.*y1.*cos(b)+(y3+2.*a).*sin(b)))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2))+y2.*sin(b).*(1./(y1.^2+y2.^2+y3.^2).*cos(b).*y1./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y1-((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1-sin(b))-1./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y1+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))))./pi./(1-nu)+1./4.*((2-2.*nu).*(y2./y1.^2./(1+y2.^2./y1.^2)-y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*cos(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+(y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*y1-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*(2.*y1.*cos(b)+(y3+2.*a).*sin(b)))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2))-(2-2.*nu).*y2.*sin(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))-(2-2.*nu).*y2.*sin(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1-y2.*(y3+a).*sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*y1-y2.*(y3+a).*sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+y2.*(y3+a).*sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1-2.*a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y1))./pi./(1-nu));
e23 = 1./2.*B1.*(1./8.*((1-2.*nu).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-1)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3)+(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-cos(b).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-cos(b))./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))))-y2.^2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).*y3-1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-1)-1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*y3+4.*a)-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)-cos(b).*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y3-1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-cos(b))-1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(2.*y3+4.*a)-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b)))))./pi./(1-nu)+1./4.*((1-2.*nu).*(((2-2.*nu).*cot(b).^2-nu).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-((2-2.*nu).*cot(b).^2+1-2.*nu).*cos(b).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))+(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(y1.*cot(b).*(1-2.*nu-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+nu.*(y3+2.*a)-a+y2.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)-(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(1./2.*a.*y1.*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)+nu-y2.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)-1./2.*y2.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a))-(1-2.*nu).*sin(b).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+(1-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+1./2.*(1-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1.*cot(b)+3./2.*a.*y1.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*(2.*y3+4.*a)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-2.*nu+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((1-2.*nu).*y1.*cot(b)-a)+y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2))-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(-2.*nu+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((1-2.*nu).*y1.*cot(b)-a)+y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((1-2.*nu).*y1.*cot(b)-a).*(2.*y3+4.*a)-1./2.*y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(2.*y3+4.*a)-y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)-1./2.*y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a.*(2.*y3+4.*a)-3./2.*a.*y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*(2.*y3+4.*a))+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b).^2-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((1-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)+a.*cos(b))+a.*(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).^2-a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)))-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b).^2-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((1-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)+a.*cos(b))+a.*(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).^2-a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a))).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((1-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)+a.*cos(b)).*(2.*y3+4.*a)-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1-2.*nu).*sin(b).*cot(b)+a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)+a.*(y3+2.*a).*sin(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)-3./2.*a.*(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*(2.*y3+4.*a)+1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).^2-a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)).*(2.*y3+4.*a)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(y2.^2.*cos(b).^2-a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(-a.*sin(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)+1./2.*a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*(2.*y3+4.*a)-a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1))))./pi./(1-nu))+1./2.*B2.*(1./8.*((2-2.*nu).*(y2./(y1.*cos(b)-y3.*sin(b)).^2.*sin(b)./(1+y2.^2./(y1.*cos(b)-y3.*sin(b)).^2)+(y2./(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).*y3+y2.*(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b).^2./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2.*y1)./(1+y2.^2.*(y1.^2+y2.^2+y3.^2).*sin(b).^2./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2)-y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*sin(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+(1./2.*y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*(2.*y3+4.*a)-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*y1)./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2))+y1.*y2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).*y3-1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-1)-1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*y3+4.*a)-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1))-y2.*(-sin(b)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-(y1.*cos(b)-y3.*sin(b))./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y3-(y1.*cos(b)-y3.*sin(b))./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-cos(b))+sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-1./2.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(2.*y3+4.*a)-(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))))./pi./(1-nu)+1./4.*((2-2.*nu).*(1-2.*nu).*(-y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*sin(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+(1./2.*y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*(2.*y3+4.*a)-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*y1)./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2)).*cot(b).^2-(1-2.*nu).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((-1+2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)+(1-2.*nu).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-1./2.*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a).*cot(b)-y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)-1./2.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a))+(1-2.*nu).*y2.*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+1./2.*(1-2.*nu).*y2.*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b).*(2.*y3+4.*a)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2.*cot(b)+3./2.*a.*y2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*(2.*y3+4.*a)+y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((1-2.*nu).*cot(b)-2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)))-1./2.*y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((1-2.*nu).*cot(b)-2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))).*(2.*y3+4.*a)-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((1-2.*nu).*cot(b)-2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)+y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)+1./2.*a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)).*(2.*y3+4.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)-1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)))+y2.*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*((-2+2.*nu).*cos(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)./cos(b))-1./2.*y2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*((-2+2.*nu).*cos(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)./cos(b)).*(2.*y3+4.*a)-y2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*((-2+2.*nu).*cos(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)./cos(b)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+y2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*((1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))-1./2.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b).*(2.*y3+4.*a)+a./(y1.^2+y2.^2+(y3+2.*a).^2)./cos(b)-a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^2./cos(b).*(2.*y3+4.*a)))./pi./(1-nu))+1./2.*B3.*(1./8.*((1-2.*nu).*sin(b).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-cos(b))./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-y2.^2.*sin(b).*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y3-1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-cos(b))-1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(2.*y3+4.*a)-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))))./pi./(1-nu)+1./4.*((1-2.*nu).*(-sin(b).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)+1./2.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)+sin(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-(y1.*cos(b)+(y3+2.*a).*sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))-1./2.*(y1.*cos(b)+(y3+2.*a).*sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a))+y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(a./(y1.^2+y2.^2+(y3+2.*a).^2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))-1./2.*y1.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(a./(y1.^2+y2.^2+(y3+2.*a).^2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)).*(2.*y3+4.*a)+y1.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-a./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*(2.*y3+4.*a)-1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1))-1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(sin(b).*(cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).*sin(b)-a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)))+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(sin(b).*(cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).*sin(b)-a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a))).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./2.*sin(b).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)+sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-1./2.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(1+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*(2.*y3+4.*a)+(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(a./(y1.^2+y2.^2+(y3+2.*a).^2)-a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*(2.*y3+4.*a))+1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).*sin(b)-a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)).*(2.*y3+4.*a)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(y2.^2.*cos(b).*sin(b)-a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(-a.*sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)+1./2.*a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*(2.*y3+4.*a)-a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1))))./pi./(1-nu))+1./2.*B1.*(1./8.*(1./(y1.^2+y2.^2+y3.^2).^(1./2)-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-cos(b).*(((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))))./pi./(1-nu)+1./8.*y2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2).*y2+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2-cos(b).*(1./(y1.^2+y2.^2+y3.^2).*cos(b).*y2./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y2-((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y2+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*y2))./pi./(1-nu)+1./4.*((2-2.*nu).*((1-2.*nu).*(-1./y1./(1+y2.^2./y1.^2)+1./(y1.*cos(b)+(y3+2.*a).*sin(b))./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))+y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))-2.*y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*cos(b))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2)).*cot(b)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-y2.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y2.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)-cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+y2.^2.*cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y2.^2.*cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2))+(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*nu./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)+a./(y1.^2+y2.^2+(y3+2.*a).^2))-y2.^2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*nu./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)+a./(y1.^2+y2.^2+(y3+2.*a).^2))+y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-2.*nu./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-2.*a./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y2)+(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1-2.*nu-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-y2.^2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1-2.*nu-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-y2.^2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1-2.*nu-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))+y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2+2.*a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y2))./pi./(1-nu))+1./2.*B2.*(1./8.*((-1+2.*nu).*sin(b).*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y2./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-y1.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2).*y2+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2)+(y1.*cos(b)-y3.*sin(b))./(y1.^2+y2.^2+y3.^2).*cos(b).*y2./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-(y1.*cos(b)-y3.*sin(b)).*((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y2-(y1.*cos(b)-y3.*sin(b)).*((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*y2-(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y2+(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*y2)./pi./(1-nu)+1./4.*((-2+2.*nu).*(1-2.*nu).*cot(b).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))+(2-2.*nu).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2+(2-2.*nu).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2-(2-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-(2-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2-(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((1-2.*nu).*cot(b)-2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2)).*y2+(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2+2.*a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y2)+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b).*sin(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(sin(b)-(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2)-(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y2.*cot(b).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))).*y2+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(sin(b)-(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2)-(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))).*y2+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y2-(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y2+(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*y2)))./pi./(1-nu))+1./2.*B3.*(1./8.*((2-2.*nu).*(1./(y1.*cos(b)-y3.*sin(b))./(1+y2.^2./(y1.*cos(b)-y3.*sin(b)).^2)+((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b))+y2.^2./(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b))-2.*y2.^2.*(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2.*cos(b))./(1+y2.^2.*(y1.^2+y2.^2+y3.^2).*sin(b).^2./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2)-1./(y1.*cos(b)+(y3+2.*a).*sin(b))./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))+y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))-2.*y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*cos(b))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2))+sin(b).*(((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))+y2.*sin(b).*(1./(y1.^2+y2.^2+y3.^2).*cos(b).*y2./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y2-((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y2+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*y2))./pi./(1-nu)+1./4.*((2-2.*nu).*(-1./y1./(1+y2.^2./y1.^2)+1./(y1.*cos(b)+(y3+2.*a).*sin(b))./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))+y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))-2.*y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*cos(b))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2))+(2-2.*nu).*sin(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-(2-2.*nu).*y2.^2.*sin(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-(2-2.*nu).*y2.^2.*sin(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)+(y3+a).*sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-y2.^2.*(y3+a).*sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-y2.^2.*(y3+a).*sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))+y2.*(y3+a).*sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2-2.*a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y2))./pi./(1-nu));
end
%====================================================
%% I_CalcStrain.m
function [A1_Ssn,A1_Sdn]=I_CalcStrain(Slip,sSsn,sSdn,dSsn,dSdn)
n=size(Slip,1);
SSlip=Slip(1:n,1:2);
slip=reshape(SSlip,2*n,1);
sdn=[sSsn dSsn; sSdn dSdn];

whos

A1_S=sdn*slip;
A1_Ssn=A1_S(1:n);
A1_Sdn=A1_S(n+1:2*n);
end

%% I_CalcStrainF.m
function [F_Ssn,F_Sdn]=I_CalcStrainF(F_Slip,sSsn,sSdn,dSsn,dSdn)

n=size(F_Slip,1);
FF_Slip=F_Slip(1:n,1:2);
slip=reshape(FF_Slip,2*n,1);
sdn=[sSsn dSsn; sSdn dSdn];
F_S=sdn*slip;
F_Ssn=F_S(1:n);
F_Sdn=F_S(n+1:2*n);
end

%====================================================
%% OutofAsperity01.m
function [A1_Slip]=OutofAsperity01(A1_Ssn,A1_Sdn,Slip,sSsn,sSdn,dSsn,dSdn)
index_out_asp=(Slip(:,1).^2+Slip(:,2).^2)==0;
%O_sSsn=sSsn(~index_asp,~index_asp);
%O_sSdn=sSdn(~index_asp,~index_asp);
%O_dSsn=dSsn(~index_asp,~index_asp);
%O_dSdn=dSdn(~index_asp,~index_asp);
%O_Ssn=A_Ssn(~index_asp);
%O_Sdn=A_Sdn(~index_asp);
%
k=sum(index_out_asp);
d=[A1_Ssn(index_out_asp);A1_Sdn(index_out_asp)];
G=[sSsn(index_out_asp,index_out_asp),dSsn(index_out_asp,index_out_asp);...
   sSdn(index_out_asp,index_out_asp),dSdn(index_out_asp,index_out_asp)];
m=G\d;
A1_Slip=Slip;
A1_Slip(index_out_asp,1:2)=-[m(1:k) m(k+1:end)];

%{
k=0;
for i=1:n
normSlip=norm(Slip(i,:));
 if normSlip==0
  k=k+1;
  k_sSsn(k,:)=sSsn(i,:);
  k_sSdn(k,:)=sSdn(i,:);
  k_dSsn(k,:)=dSsn(i,:);
  k_dSdn(k,:)=dSdn(i,:);
  O_Ssn(k,1)=A_Ssn(i);
  O_Sdn(k,1)=A_Sdn(i);
 else
 end
end

m=0;

for j=1:n
normSlip=norm(Slip(j,:));
 if normSlip==0
  m=m+1;
  O_sSsn(:,m)=k_sSsn(:,j);
  O_sSdn(:,m)=k_sSdn(:,j);
  O_dSsn(:,m)=k_dSsn(:,j);
  O_dSdn(:,m)=k_dSdn(:,j);
 else
 end
end
%}
end
%% OutofAsperity01.m
function [FO_Ssn,FO_Sdn]=OutofAsperity11(F_Ssn,F_Sdn,Slip)
n=length(Slip);
k=0;
for i=1:n
normSlip=norm(Slip(i,:));
 if normSlip==0
  k=k+1;
 
  FO_Ssn(k,1)=F_Ssn(i);
  FO_Sdn(k,1)=F_Sdn(i);
 else
 end
end
end

%====================================================
%% SasaEular2Velo.m
function [vnorm,v,NARFA]=SasaEular2Velo(ll)

Fid=fopen('/home_tmp/sasajima/DATA/Pole_data/PAC-OKH_MORVEL56.dat','r');
Epole=textscan(Fid,'%f %f %f');
fclose(Fid);
Epole=cell2mat(Epole);

%Fid2=fopen('/home_tmp/sasajima/DATA/boso/IA-KA_vll.txt','w');

R(1,1)=6371000;%[m]
%a=6378137.0;
%b=6356752.3;
%aabb=(a.^2)./(b.^2);

diglonP(1,1)=Epole(1,1);%longutitude of Eular pole [digree]%
diglatP(1,1)=Epole(1,2);%latitude of Eular pole [digree]%
digwMa(1,1)=Epole(1,3);%Rotation velocity [digree/Ma]%
radlonP=diglonP./180.*pi;
radlatP(1,1)=diglatP(1,1)./180.*pi;
radwMa=digwMa./180.*pi;

radlonX(:,1)=ll(:,1)./180.*pi;%lon of observed point[radian]
radlatX(:,1)=ll(:,2)./180.*pi;%lat of observed point[radian]
diglonX(:,1)=ll(:,2);%lon of observed poiint [digree]

NARFA(:,1)=radlonX(:,1);
NARFA(:,2)=radlatX(:,1);

CE(1,1)=0.5.*pi-radlatP(1,1);%[radian],scolar
CB(:,1)=0.5.*pi-radlatX(:,1);%[radian],vector
RMD(:,1)=radlonP-radlonX(:,1);%[radian],vector

cCE(1,1)=cos(CE);
sCE(1,1)=sin(CE);

n=length(ll);
v=zeros(n,2);
vnorm=zeros(n,1);

for i=1:n
    cCB(1,1)=cos(CB(i,1));
    sCB(1,1)=sin(CB(i,1));
    cRMD(1,1)=cos(RMD(i,1));
    
    cDLT(1,1)=cCB(1,1).*cCE(1,1)+sCB(1,1).*sCE(1,1).*cRMD(1,1);
    DLT(1,1)=acos(cDLT(1,1));
    
    sDLT(1,1)=sin(DLT(1,1));
    L(1,1)=sDLT(1,1).*R(1,1);%distance between obs point and Eular Vector [m]
    
    vmMa(1,1)=L(1,1).*radwMa(1,1);%velocity[m/Ma]
    vnorm(i,1)=vmMa(1,1)./10000;%velocity[cm/year]
    
    cARFA(1,1)=(cCE(1,1)-(cCB(1,1).*cDLT(1,1)))./(sCB(1,1).*sDLT(1,1));
    ARFA(1,1)=acos(cARFA(1,1));%[radian]
    
    dltdiglon(1,1)=diglonP(1,1)-diglonX(i,1);%[digree]
    
    if dltdiglon>0
        NARFA(i,3)=0.5.*pi-ARFA(1,1);%countour? clock wise from North [radian]
    else
        NARFA(i,3)=0.5.*pi+ARFA(1,1);%countour? clock wise from North[radian]
    end
    
    v(i,1)=vnorm(i,1).*cos(NARFA(i,3));%vlat;N=plus[cm/year]
    v(i,2)=vnorm(i,1).*sin(NARFA(i,3));%vlon;E=plus[cm/year]
    
    %digNARFA(i,1)=NARFA(i,3).*180./pi;%clock wise from North[digree]
    
    %fprintf(Fid2,'%15.8f %15.8f %15.8f %15.8f\n',v(i,1),v(i,2),vnorm(i,1),digNARFA(i,1));
    
end

%fclose(Fid2);
end

%% Sasa_ll2xyz.m
function [trixyzC,trixyz3,sxyz]=trill2trixyz(triC,tri3,sll,ALAT0,ALON0)
% Convert triC, tri3, and sll from LonLat to local XY cordinate.
% Original coded by Sasajima.
% Modified by H. Kimura in 2018/11/13.

nn=size(triC,1);
% mm=size( sll,1);
% trixyz3=zeros(nn,3,3);
% trixyzC=zeros(nn,3);

%*** Convert of sll
[sxyz(:,2),sxyz(:,1)]=PLTXY(sll(:,1),sll(:,2),ALAT0,ALON0);

%*** Convert of tri3
tri3lon=[tri3(:,1);tri3(:,4);tri3(:,7)];
tri3lat=[tri3(:,2);tri3(:,5);tri3(:,8)];
tri3dep=[tri3(:,3);tri3(:,6);tri3(:,9)];
[tmpxyz3(:,2),tmpxyz3(:,1)]=PLTXY(tri3lat,tri3lon,ALAT0,ALON0);
tmpxyz3(:,3)=1e3.*tri3dep;
trixyz3=[tmpxyz3(1:nn,:) tmpxyz3(1+nn:2*nn,:) tmpxyz3(1+2*nn:3*nn,:)];

%*** Convert of triC
[trixyzC(:,2),trixyzC(:,1)]=PLTXY(triC(:,2),triC(:,1),ALAT0,ALON0);
trixyzC(:,3)=1e3.*triC(:,3);

end
%====================================================

%% FSlip2xyll.m
function [xyFSlip]=FSlip2xyll(F_Slip,sitaS)
xyFSlip(:,1)=F_Slip(:,1).*cos(sitaS(:))+F_Slip(:,2).*sin(sitaS(:));
xyFSlip(:,2)=-F_Slip(:,1).*sin(sitaS(:))+F_Slip(:,2).*cos(sitaS(:));
end
%====================================================
%% defineslipQ.m
function [Slip]=defineSlipQ(triC,sitaS)
% A1=Rectangle asperity_1%
% p=minimum X and max Y, q=max X and minimum Y%
% y=+ wa south%

% 2011 Tohoku-oki
A1SV=1.8810;
A1plon=143;
A1plat=41;
A1qlon=145;
A1qlat=38;
A1Slip=0.090;
% %Hokkaido 500year
% A1SV=1.8810;  % unknown
% A1plon=145.0;
% A1plat=42.2;
% A1qlon=146.5;
% A1qlat=41.25;
% A1Slip=9.2;
% %1896 Meiji Sanriku-oki
% A1SV=1.8810;  % unknown
% A1plon=143.7;
% A1plat=39.6;
% A1qlon=144.15;
% A1qlat=38.8;
% A1Slip=8.85;
% %2003 Tokachi-oki
% A1SV=1.8810;  % unknown
% A1plon=143.65;
% A1plat=42.2;
% A1qlon=144.2;
% A1qlat=41.7;
% A1Slip=8.5;
A1xSlip=-A1Slip.*sin(A1SV);
A1ySlip=A1Slip.*cos(A1SV);

n=size(triC,1);
Slip=zeros(n,2);
define=zeros(n,1);
for i=1:n
  if triC(i,1)<A1plon,Slip(i,:)=[0,0];
  elseif triC(i,1)>A1qlon,Slip(i,:)=[0,0];
  elseif triC(i,2)<A1qlat,Slip(i,:)=[0,0];
  elseif triC(i,2)>A1plat,Slip(i,:)=[0,0];
  else
    A1sSlip=cos(sitaS(i)).*A1xSlip-sin(sitaS(i)).*A1ySlip;
    A1dSlip=sin(sitaS(i)).*A1xSlip+cos(sitaS(i)).*A1ySlip;
    Slip(i,:)=[A1sSlip,A1dSlip];
    define(i,1)=1;
  end
end

end
%====================================================
%% Make patches of locked asperity
function [Slip]=MakeSlipPatch(triC,sitaS)
% Read file of locked patch definition
fid=fopen(prm.patch,'r');
np = 0;
nf = 0;
while 1
  tline = fgetl(fid);
  if ~ischar(tline)
    break
  end
  nf  = nf + 1;
  tmp = strsplit(tline);
  tmp = str2double(tmp);
  if isnan(tmp)
    np = np + 1;
    nf = 0;
  else
    Slip.patch(np+1).lon(nf) = tmp(1);
    Slip.patch(np+1).lat(nf) = tmp(2);
  end
end
fclose(fid);

% Generate locked patch


% A1=Rectangle asperity_1%
% p=minimum X and max Y, q=max X and minimum Y%
% y=+ wa south%

% 2011 Tohoku-oki
A1SV=1.8810;
A1plon=143;
A1plat=41;
A1qlon=145;
A1qlat=38;
A1Slip=0.090;
A1xSlip=-A1Slip.*sin(A1SV);
A1ySlip=A1Slip.*cos(A1SV);

n=size(triC,1);
Slip=zeros(n,2);
define=zeros(n,1);
for i=1:n
  if triC(i,1)<A1plon,Slip(i,:)=[0,0];
  elseif triC(i,1)>A1qlon,Slip(i,:)=[0,0];
  elseif triC(i,2)<A1qlat,Slip(i,:)=[0,0];
  elseif triC(i,2)>A1plat,Slip(i,:)=[0,0];
  else
    A1sSlip=cos(sitaS(i)).*A1xSlip-sin(sitaS(i)).*A1ySlip;
    A1dSlip=sin(sitaS(i)).*A1xSlip+cos(sitaS(i)).*A1ySlip;
    Slip(i,:)=[A1sSlip,A1dSlip];
    define(i,1)=1;
  end
end

end

%% Make Green's function of halfspace elastic strain for triangular meshes
function [tri] = MakeGreenFunction(blk,obs)
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
tri(1).cf    = ones(3*blk(1).nb,1);
tri(1).inv   = zeros(3*blk(1).nb,1);
for nb1 = 1:blk(1).nblock
  [blk(nb1).localx,blk(nb1).localy] = PLTXY(blk(nb1).lat,blk(nb1).lon,alat,alon);
  for nb2 = nb1+1:blk(1).nblock
    nf = size(blk(1).bound(nb1,nb2).blon,1);
    tri(1).bound(nb1,nb2).clat = [];
    tri(1).bound(nb1,nb2).clon = [];
    tri(1).bound(nb1,nb2).cdep = [];
    if nf ~= 0
      tri(1).bound(nb1,nb2).gstr = zeros(3*nd,nf);
      tri(1).bound(nb1,nb2).gdip = zeros(3*nd,nf);
      tri(1).bound(nb1,nb2).gtns = zeros(3*nd,nf);
%
      fprintf('==================\n Block %2i : Block %2i \n Number of TRI sub-faults : %4i \n',nb1,nb2,nf)
%
      triclon = mean(blk(1).bound(nb1,nb2).blon,2);
      triclat = mean(blk(1).bound(nb1,nb2).blat,2);
      tricdep = mean(blk(1).bound(nb1,nb2).bdep,2);
      [tricx,tricy] = PLTXY(triclat,triclat,alat,alon);
      tricz = -tricdep;

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
        tri(1).bound(nb1,nb2).phi(n)    = phi;
        tri(1).bound(nb1,nb2).theta(n)  = theta;
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
      tri(1).nb = tri(1).nb+nf;
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
    tri(1).cf(3*tri(1).nb+nf+n) = 1/sqrt(dp(1)^2+dp(2)^2);  % 1=sqrt(dp(1)^2+dp(2)^2+dp(3)^2): norm of dp
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
    tri(1).inv(3*tri(1).nb     +n) =  1;
    tri(1).inv(3*tri(1).nb+  nf+n) =  1;
    tri(1).inv(3*tri(1).nb+2*nf+n) =  0;
  case 2
    tri(1).inv(3*tri(1).nb     +n) = -1;
    tri(1).inv(3*tri(1).nb+  nf+n) = -1;
    tri(1).inv(3*tri(1).nb+2*nf+n) =  0;
  case 0
    trixc   = mean(trix);
    triyc   = mean(triy);
    [in,on] = inpolygon(trixc,triyc,blk(nb1).localx,blk(nb1).localy);
    if in==1 && on~=1
      tri(1).inv(3*tri(1).nb     +n) =  1;
      tri(1).inv(3*tri(1).nb+  nf+n) =  0;
      tri(1).inv(3*tri(1).nb+2*nf+n) =  1;
    elseif in==1 && on==1
      tric = [trixc triyc 0];
      uv   = [0 0 1];
      nv   = cross(uv,tri(1).bound(nb1,nb2).st(n,:));
      cnv  = tric+nv;
      if inpolygon(cnv(1),cnv(2),blk(nb1).localx,blk(nb1).localy)==1
        tri(1).inv(3*tri(1).nb     +n) =  1;
        tri(1).inv(3*tri(1).nb+  nf+n) =  0;
        tri(1).inv(3*tri(1).nb+2*nf+n) =  1;
      else
        tri(1).inv(3*tri(1).nb     +n) = -1;
        tri(1).inv(3*tri(1).nb+  nf+n) =  0;
        tri(1).inv(3*tri(1).nb+2*nf+n) = -1;
      end
    else
      tri(1).inv(3*tri(1).nb     +n) =  -1;
      tri(1).inv(3*tri(1).nb+  nf+n) =   0;
      tri(1).inv(3*tri(1).nb+2*nf+n) =  -1;
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

%% make Green's function of displacement
function [Gu]=makeGreenDisp(obs,trixyz3,tri)
% This function make Green's function of displacement (Gu).
% Input
%  xyz     : site location.
%  trixyz3 : trimesh coordianate.
% Output
%  Gu.st   : surface displacement due to strike slip on a fault.
%  Gu.dp   : surface displacement due to dip slip on a fault.
%  Gu.ts   : surface displacement due to tensile slip on a fault.
%  *Each matrix (Gu.*) contain 3 x NOBS by NFLT elements.

pr  =0.25; % Poisson's ratio
nobs=size(obs.x,2);
% nobs=size(xyz,1);
nflt=size(trixyz3,1);
sx  =obs.x';
sy  =obs.y';
sz  =obs.z';
% sx  =xyz(:,1);
% sy  =xyz(:,2);
% sz  =zeros(nobs,1);

for n=1:nflt
    trix=[trixyz3(n,1),trixyz3(n,4),trixyz3(n,7)];
    triy=[trixyz3(n,2),trixyz3(n,5),trixyz3(n,8)];
    triz=[trixyz3(n,3),trixyz3(n,6),trixyz3(n,9)];
%     [U] = CalcTriDisps(sx, sy, sz, x, y, z, pr, ss, ts, ds);
    [U] = CalcTriDisps(sx, sy, sz, trix, triy, triz, pr, 1, 0, 0); % Strike
    Gu.st(1:3:3*nobs,n)= U.x;
    Gu.st(2:3:3*nobs,n)= U.y;
    Gu.st(3:3:3*nobs,n)=-U.z;
    [U] = CalcTriDisps(sx, sy, sz, trix, triy, triz, pr, 0, 1, 0); % Tensile
    Gu.ts(1:3:3*nobs,n)= U.x;
    Gu.ts(2:3:3*nobs,n)= U.y;
    Gu.ts(3:3:3*nobs,n)=-U.z;
    [U] = CalcTriDisps(sx, sy, sz, trix, triy, triz, pr, 0, 0, 1); % Dip
    Gu.dp(1:3:3*nobs,n)= U.x;
    Gu.dp(2:3:3*nobs,n)= U.y;
    Gu.dp(3:3:3*nobs,n)=-U.z;
end

end
%% make Green's function of strain
function [Gs]=makeGreenStrain(trixyzC,trixyz3,sitaS,sitaD,normVec)
% This function make Green's function of elastic strain (Gs).
% Input
%  trixyzC : center of trimesh.
%  trixyz3 : trimesh coordianate.
%  sitaS   : angle of fault strike from X-axis.
%  sitaD   : angle of fault dip from XY-plane.
%  normVec : normalized normal vectors of faults.
% Output
%  Gs.st   : strain due to strike slip on a fault.
%  Gs.dp   : strain due to dip slip on a fault.
%  Gs.ts   : strain due to tensile slip on a fault.
%  *Each matrix (Gs.*) contain 6 x NFLT by NFLT elements.
pr  =0.25; % Poisson's ratio
nflt=size(trixyz3,1);
sx=trixyzC(:,1)+normVec(:,1);
sy=trixyzC(:,2)+normVec(:,2);
sz=trixyzC(:,3)+normVec(:,3);

for n=1:nflt
    trix=[trixyz3(n,1),trixyz3(n,4),trixyz3(n,7)];
    triy=[trixyz3(n,2),trixyz3(n,5),trixyz3(n,8)];
    triz=[trixyz3(n,3),trixyz3(n,6),trixyz3(n,9)];
%     [S] = CalcTriStrains(sx, sy, sz, x, y, z, pr, ss, ts, ds);
    [S] = CalcTriStrains(sx, sy, sz, trix, triy, triz, pr, 1, 0, 0); % Strike
    Gs.st(1:6:6*nflt,n)=S.xx;
    Gs.st(2:6:6*nflt,n)=S.xy;
    Gs.st(3:6:6*nflt,n)=S.xz;
    Gs.st(4:6:6*nflt,n)=S.yy;
    Gs.st(5:6:6*nflt,n)=S.yz;
    Gs.st(6:6:6*nflt,n)=S.zz;
    [S] = CalcTriStrains(sx, sy, sz, trix, triy, triz, pr, 0, 1, 0); % Tensile
    Gs.ts(1:6:6*nflt,n)=S.xx;
    Gs.ts(2:6:6*nflt,n)=S.xy;
    Gs.ts(3:6:6*nflt,n)=S.xz;
    Gs.ts(4:6:6*nflt,n)=S.yy;
    Gs.ts(5:6:6*nflt,n)=S.yz;
    Gs.ts(6:6:6*nflt,n)=S.zz;
    [S] = CalcTriStrains(sx, sy, sz, trix, triy, triz, pr, 0, 0, 1); % Dip
    Gs.dp(1:6:6*nflt,n)=S.xx;
    Gs.dp(2:6:6*nflt,n)=S.xy;
    Gs.dp(3:6:6*nflt,n)=S.xz;
    Gs.dp(4:6:6*nflt,n)=S.yy;
    Gs.dp(5:6:6*nflt,n)=S.yz;
    Gs.dp(6:6:6*nflt,n)=S.zz;
end
[Gs]=trans_xyz2strdip_old(Gs,sitaS,sitaD);
end

%%
function [tri] = trans_xyz2strdip(tri,nb1,nb2,n,nf)
% This function transforms a strain tensor from xyz to fault strike-dip.
% Strain of strike direction on the fault corresponds to ezx',
% strain of dip direction on the fault corresponds to ezy', and
% strain of tensile direction on the fault corresponds to ezz'.

c1=cos(tri(1).bound(nb1,nb2).phi(n));
s1=sin(tri(1).bound(nb1,nb2).phi(n));
c2=cos(tri(1).bound(nb1,nb2).theta(n));
s2=sin(tri(1).bound(nb1,nb2).theta(n));

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

%% Transform tensor from xyz to strike-dip
function [Gs]=trans_xyz2strdip_old(Gs,sitaS,sitaD)
% This function transforms a strain tensor from xyz to fault strike-dip.
% 
% Output
% Gs.stst : strain of strike direction on the fault due to strike slip.
% Gs.stdp : strain of dip direction on the fault due to strike slip.
% Gs.stts : strain of tensile direction on the fault due to strike slip.
% Gs.dpst : strain of strike direction on the fault due to dip slip.
% Gs.dpdp : strain of dip direction on the fault due to dip slip.
% Gs.dpts : strain of tensile direction on the fault due to dip slip.
% Gs.tsst : strain of strike direction on the fault due to tensile slip.
% Gs.tsdp : strain of dip direction on the fault due to tensile slip.
% Gs.tsts : strain of tensile direction on the fault due to tensile slip.
% 
% Note that strain of strike direction on the fault corresponds to ezx',
% strain of dip direction on the fault corresponds to ezy', and
% strain of tensile direction on the fault corresponds to ezz'.

% c1=cos(sitaS);
% s1=sin(sitaS);
% c2=cos(sitaD);
% s2=sin(sitaD);
c1=repmat(cos(sitaS)',1,size(sitaS,2));
s1=repmat(sin(sitaS)',1,size(sitaS,2));
c2=repmat(cos(sitaD)',1,size(sitaD,2));
s2=repmat(sin(sitaD)',1,size(sitaD,2));

[sst,sdp,sts]=calctrans(Gs.st,c1,s1,c2,s2); % response to strike slip
Gs.stst=sst;
Gs.stdp=sdp;
Gs.stts=sts;
[sst,sdp,sts]=calctrans(Gs.ts,c1,s1,c2,s2); % response to tensile slip
Gs.tsst=sst;
Gs.tsdp=sdp;
Gs.tsts=sts;
[sst,sdp,sts]=calctrans(Gs.dp,c1,s1,c2,s2); % response to dip slip
Gs.dpst=sst;
Gs.dpdp=sdp;
Gs.dpts=sts;

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
%% make_test_trill.m
%====================================================
function [triC,tri3,tri,sll]=make_test_trill
%by Ryohei Sasajima 2013/12/23
%====================================================
%====================================================
Fid0=fopen('/home/sasajima/Dropbox/yellow/PACdepth201312.txt','r');
dep_main=textscan(Fid0,'%f %f %f');
fclose(Fid0);
dep_main=cell2mat(dep_main);
%====================================================
Fid00=fopen('/home/sasajima/Dropbox/yellow/depth0_1402.txt','r');
dep_sub=textscan(Fid00,'%f %f %f');
fclose(Fid00);
dep_sub=cell2mat(dep_sub);
%===================================================
%====================================================
Fid1=fopen('/home/sasajima/Dropbox/yellow/PAC_blue.txt','r');
bound=textscan(Fid1,'%f %f');
fclose(Fid1);
bound=cell2mat(bound);

Fid2=fopen('/home/sasajima/Dropbox/yellow/PAC_blue.txt','r');
PACblue=textscan(Fid2,'%f %f');
fclose(Fid2);
PACblue=cell2mat(PACblue);

Fid3=fopen('/home/sasajima/Dropbox/yellow/PAC_green.txt','r');
PACgreen=textscan(Fid3,'%f %f');
fclose(Fid3);
PACgreen=cell2mat(PACgreen);

Fid4=fopen('/home_tmp/sasajima/DATA/blue.dat','r');
mblue=textscan(Fid4,'%f %f');
fclose(Fid4);
mblue=cell2mat(mblue);
bb=size(mblue,1);

Fid5=fopen('/home_tmp/sasajima/DATA/green.dat','r');
mgreen=textscan(Fid5,'%f %f');
fclose(Fid5);
mgreen=cell2mat(mgreen);
gg=size(mgreen,1);

Fid7=fopen('/home/sasajima/Dropbox/yellow/PAC_red.txt','r');
PACred=textscan(Fid7,'%f %f');
fclose(Fid7);
PACred=cell2mat(PACred);

Fid8=fopen('/home/sasajima/Dropbox/yellow/PAC_black.txt','r');
PACblack=textscan(Fid8,'%f %f');
fclose(Fid8);
PACblack=cell2mat(PACblack);

Fid9=fopen('/home_tmp/sasajima/DATA/red.dat','r');
mred=textscan(Fid9,'%f %f');
fclose(Fid9);
mred=cell2mat(mred);
rr=size(mred,1);

Fid10=fopen('/home_tmp/sasajima/DATA/black.dat','r');
mblack=textscan(Fid10,'%f %f');
fclose(Fid10);
mblack=cell2mat(mblack);
blbl=size(mblack,1);

Fid11=fopen('/home_tmp/sasajima/DATA/m_all.dat','r');
mall=textscan(Fid11,'%f %f');
fclose(Fid11);
mall=cell2mat(mall);
alal=size(mall,1);


s=0;
for bl=1:blbl
    Blue=inpolygon(mblack(bl,1),mblack(bl,2),PACblue(:,1),PACblue(:,2));
    if Blue==0
    else
        s=s+1;
        slon(s,1)=mblack(bl,1);
        slat(s,1)=mblack(bl,2);
    end
end

%{

 for b=1:bb;

 Blue=inpolygon(mblue(b,1),mblue(b,2),PACblue(:,1),PACblue(:,2));
 BRed=inpolygon(mblue(b,1),mblue(b,2),PACred(:,1),PACred(:,2));

  if (Blue==1)&&(BRed==0)
   s=s+1;
   slon(s,1)=mblue(b,1);
   slat(s,1)=mblue(b,2);
  else
  end
 end

 sblue=s-sred

 for g=1:gg;

  Green=inpolygon(mgreen(g,1),mgreen(g,2),PACgreen(:,1),PACgreen(:,2));
 GBlue=inpolygon(mgreen(g,1),mgreen(g,2),PACblue(:,1),PACblue(:,2));

  if (Green==1)&&(GBlue==0)
   s=s+1;
   slon(s,1)=mgreen(g,1);
   slat(s,1)=mgreen(g,2);
  else
  end
 end

 sgreen=s-sred-sblue

 for bl=1:blbl;

  Black=inpolygon(mblack(bl,1),mblack(bl,2),PACblack(:,1),PACblack(:,2));
 BlGreen=inpolygon(mblack(bl,1),mnlack(bl,2),PACgreen(:,1),PACgreen(:,2));

  if (Black==1)&&(BlGreen==0)
   s=s+1;
   slon(s,1)=mblack(bl,1);
   slat(s,1)=mblack(bl,2);
  else
  end
 end

 sblack=s-sred-sblue-sgreen


 for al=1:alal;

  All=inpolygon(mall(al,1),mall(al,2),bound(:,1),bound(:,2));
 AlBlack=inpolygon(mall(al,1),mall(al,2),PACblack(:,1),PACblack(:,2));

  if (All==1)&&(AlBlack==0)
   s=s+1;
   slon(s,1)=mall(al,1);
   slat(s,1)=mall(al,2);
  else
  end
 end

 soutblack=s-sred-sblue-sgreen-sblack
 sall=s
 pause

%}

sll(:,1:2)=[slon(:,1),slat(:,1)];

%ll=[s.lon,s.lat];
%plot(s.lon,s.lat,'.')
%====================================================
tri = delaunay(slon,slat);
%====================================================
ntri=length(tri);
E=scatteredInterpolant(dep_main(:,1),dep_main(:,2),dep_main(:,3),'natural');
F=scatteredInterpolant(dep_sub(:,1),dep_sub(:,2),dep_sub(:,3),'natural');

tri3=zeros(ntri,3,3);
triC=zeros(ntri,3);
nn=0;
for n=1:ntri
  lon1=slon(tri(n,1));
  lat1=slat(tri(n,1));
  dep1=-F(lon1,lat1)-E(lon1,lat1);
  lon2=slon(tri(n,2));
  lat2=slat(tri(n,2));
  dep2=-F(lon2,lat2)-E(lon2,lat2);
  lon3=slon(tri(n,3));
  lat3=slat(tri(n,3));
  dep3=-F(lon3,lat3)-E(lon3,lat3);
  
  glon=(lon1+lon2+lon3)/3;
  glat=(lat1+lat2+lat3)/3;
  gdep=(dep1+dep2+dep3)/3;
  
  IN=inpolygon(glon,glat,bound(:,1),bound(:,2));
  
  if gdep==0
  elseif IN==0
  else
    nn=nn+1;
    triC(nn,1:3)=[glon,glat,gdep];
    tri3(nn,1)=lon1;
    tri3(nn,2)=lat1;
    tri3(nn,3)=dep1;
    tri3(nn,4)=lon2;
    tri3(nn,5)=lat2;
    tri3(nn,6)=dep2;
    tri3(nn,7)=lon3;
    tri3(nn,8)=lat3;
    tri3(nn,9)=dep3;
  end
end
tri3(nn+1:ntri,:)=[];
triC(nn+1:ntri,:)=[];

Fid6=fopen('/home_tmp/sasajima/DATA/PAC_tri3.txt','w');
for w=1:nn
  for u=1:3
    fprintf(Fid6,'%10.4f %9.4f %10.4f\n',tri3(w,1,u),tri3(w,2,u),tri3(w,3,u));
  end
end
fclose(Fid6);

end
%% makexyz.m
function  [xyz]=makexyz
%by Ryohei Sasajima
%last 2014/02/02
%====================================================
%{
Fid1=fopen('/home_tmp/sasajima/DATA/GEODETIC_DATA/coordinates_F3/sitelocate.txt','r');
geonet=textscan(Fid1,'%f %f %f');
fclose(Fid1);
geonet=cell2mat(geonet);
%}

Fid2=fopen('/home_tmp/sasajima/DATA/xyz01.dat','r');
xyz01=textscan(Fid2,'%f %f');
fclose(Fid2);
xyz=cell2mat(xyz01);

end
%====================================================
%% Slip2xyll.m
function [xySlip]=Slip2xyll(Slip,sitaS)
xySlip(:,1)=Slip(:,1).*cos(sitaS(:))+Slip(:,2).*sin(sitaS(:));
xySlip(:,2)=-Slip(:,1).*sin(sitaS(:))+Slip(:,2).*cos(sitaS(:));
end
%====================================================

%% strike_dip.m
function[sitaS,sitaD,normVec]=strike_dip(trixyzC,trixyz3)
%Define strike and dip of fault%
n=size(trixyzC,1);
ca=zeros(n,3);
cb=zeros(n,3);
InormVec=zeros(n,3);
normVec=zeros(n,3);
strikeVec=zeros(n,3);
dipVec=zeros(n,3);
sitaS=zeros(n,1);
sitaD=zeros(n,1);

for i=1:n
    ca(i,:)=[trixyz3(i,1)-trixyz3(i,7),trixyz3(i,2)-trixyz3(i,8),trixyz3(i,3)-trixyz3(i,9)];
    cb(i,:)=[trixyz3(i,1)-trixyz3(i,4),trixyz3(i,2)-trixyz3(i,5),trixyz3(i,3)-trixyz3(i,6)];
    InormVec(i,:)              = cross(ca(i,:),cb(i,:),2);% direct to deeper %
    normVec(i,:)              = InormVec(i,:)./(sqrt((InormVec(i,1)).^2+(InormVec(i,2)).^2+(InormVec(i,3).^2)));
    if (normVec(i,3) < 0) % Enforce clockwise circulation
        normVec(i,:)               = -normVec(i,:);
    end
    strikeVec(i,:)              = [-sin(atan2(normVec(i,2),normVec(i,1))),cos(atan2(normVec(i,2),normVec(i,1))),0];%direct to left hand of who toward dip direction
    dipVec(i,:)                 = cross(normVec(i,:), strikeVec(i,:),2);
    %if normVec(1)==0 and normVec(2)==0
    % strikeVec=[1 0 0];
    % dipVec   =[0 1 0];
    % normVec  =[0 0 1];
    %end
    sitaS(i)=atan2(strikeVec(i,2),strikeVec(i,1)); %from x-axis to y-axis rotation [rad]%
    sitaD(i)=asin(dipVec(i,3));
end
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
