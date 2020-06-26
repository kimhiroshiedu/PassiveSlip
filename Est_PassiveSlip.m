function Est_PassiveSlip(varargin)
close all
% Estimates Euler vectors, internal strains, and asperity lines.
% Or calculate Passive back-slip rates for given asperities.
% 
% Est_PassiveSlip(estimation,method)
% Est_PassiveSlip(method)
% where estimation: 'mcmc'(default) or 'fwd'
%       method: 'mh'(default), 're'
% e.g.) Est_PassiveSlip('fwd')
%       Est_PassiveSlip('mcmc','re')
%       Est_PassiveSlip('re')

% Coded by Ryohei Sasajima final 2013/12/23
% Combined by Hiroshi Kimura     2018/11/12
% Revised by Hiroshi Kimura      2019/02/16
% Revised by Hiroshi Kimura      2019/07/09
% Revised by Hiroshi Kimura      2020/03/09
%--- 
prm.input      = 'PARAMETER/parameter_inversion_swjp.txt';
prm.optfile    = 'PARAMETER/opt_bound_par.txt';
%--
% Read Parameters.
[prm]     = ReadParameters(prm);
% Simulation mode select
[prm]     = DetermineCalcMethod(prm,varargin);
clear varargin
% Read observation data. 
[obs]     = ReadObs(prm);
% Read blocks.
[blk,obs] = ReadBlockBound(prm,obs);
% Read dipping boundaries.
[blk,prm] = ReadBoundaryTypes(blk,prm);
% Read or generate block interfaces.
[blk]     = ReadBlockInterface(blk,prm,obs);
% Show block bound map
ShowBlockBound(blk)
% Read Euler pole vectors.
[eul,prm] = ReadEulerPoles(blk,prm);
% Read internal strain flag.
[blk,prm] = ReadInternalStrain(blk,obs,prm);
% Read UDline settings.
[blk]     = ReadUDlineSettings(blk,prm);
% Read random walk line definitions.
[blk,asp] = DefRandomWalkLine(blk,prm,obs);
% Estimate rigid motion and calculate AIC.
[blk,obs] = CalcAIC(blk,obs,eul,prm);
% Generate Green's functions.
[blk,tri] = GreenTri(blk,obs,prm);
% Arrange Green's functions.
[d,G]     = GreenFunctionMatrix(blk,obs,tri);

if strcmpi(prm.method,'Forward')
  % Read asperity areas.
  [blk]   = ReadLockedPatch(blk,prm);
  % Define locking patches.
  [d]     = InitialLockingPatch(blk,tri,d);
  % Calculate pasive slip and response to surface.
  [cal]   = CalcPassiveSlip(blk,asp,tri,prm,obs,eul,d,G);
elseif strcmpi(prm.method,'MCMC.MH')
  % Metropolis-Hasting
  [cal]   = Proceed_MCMC_MH(blk,asp,tri,prm,obs,eul,d,G);
elseif strcmpi(prm.method,'MCMC.RE')
  % Replica exchange Monte Carlo
  [cal]   = Proceed_MCMC_RE(blk,asp,tri,prm,obs,eul,d,G);
end

% Save data.
SaveData(prm,blk,obs,tri,d,G,cal)

end

%% Determine simulation method
function [prm] = DetermineCalcMethod(prm,varin)
correct = 0;
while correct ~= 1
  if ~isempty(varin)
    if strcmpi(char(varin{1}),'fwd')
      prm.method = 'Forward';
      correct = 1;
    elseif strcmpi(char(varin{1}),'mcmc')
      if size(varin,2) == 1
        prm.method = 'MCMC.MH';
        correct = 1;
      elseif size(varin,2) == 2
        if sum(strcmpi(char(varin{2}),{'mh','re'})) >= 1
          prm.method = ['MCMC.',char(upper(varin{2}))];
          correct = 1;
        else
          varin{2} = input("Invalid option. Enter 'mh' or 're' > ",'s');
        end
      else
        fprintf('Too much variables. Enter options again.\n')
        varin{1} = input('option1 > ','s');
        if strcmpi(char(varin{1}),'fwd')
          prm.method = 'Forward';
          correct = 1;
        elseif strcmpi(char(varin{1}),'mcmc')
          varin{2} = input('option2 > ','s');
        elseif sum(strcmpi(char(varin{1}),{'mh','re'})) >= 1
          prm.method = ['MCMC.',char(upper(varin{1}))];
          correct = 1;
        end
      end
    elseif sum(strcmpi(char(varin{1}),{'mh','re'})) >= 1
      prm.method = ['MCMC.',char(upper(varin{1}))];
      correct = 1;
    else
      fprintf('Invalid valiable. Enter correct options.\n')
      varin{1} = input('option1 > ','s');
      if strcmpi(char(varin{1}),'fwd')
        prm.method = 'Forward';
        correct = 1;
      elseif strcmpi(char(varin{1}),'mcmc')
        varin{2} = input('option2 > ','s');
      elseif sum(strcmpi(char(varin{1}),{'mh','re'})) >= 1
        prm.method = ['MCMC.',char(upper(varin{1}))];
        correct = 1;
      end
    end
  else
    prm.method = 'MCMC.MH';
    correct = 1;
  end
end

if sum(strcmpi(prm.method,{'Forward','MCMC.MH'})) >= 1
  prm.nrep = 1;
end

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
filepole           = fscanf(fid,'%s \n',[1,1]); [~] = fgetl(fid);
filebound          = fscanf(fid,'%s \n',[1,1]); [~] = fgetl(fid);
fileinternal       = fscanf(fid,'%s \n',[1,1]); [~] = fgetl(fid);
fileudconfig       = fscanf(fid,'%s \n',[1,1]); [~] = fgetl(fid);
dirresult          = fscanf(fid,'%s \n',[1,1]); [~] = fgetl(fid);
prm.home_d = pwd;
prm.fileobs            = fullfile(prm.home_d,fileobs);
prm.dirblock           = fullfile(prm.home_d,dirblock);
prm.dirblock_interface = fullfile(prm.home_d,dirblock_interface);
prm.filepole           = fullfile(prm.home_d,filepole);
prm.filebound          = fullfile(prm.home_d,filebound);
prm.fileinternal       = fullfile(prm.home_d,fileinternal);
prm.fileUDconfig       = fullfile(prm.home_d,fileudconfig);
prm.dirresult          = fullfile(prm.home_d,dirresult);
%
prm.gpu = fscanf(fid,'%d \n',[1,1]); [~] = fgetl(fid);
prm.itr = fscanf(fid,'%d \n',[1,1]); [~] = fgetl(fid);
prm.thr = fscanf(fid,'%d \n',[1,1]); [~] = fgetl(fid);
prm.cha = fscanf(fid,'%d \n',[1,1]); [~] = fgetl(fid);
prm.kep = fscanf(fid,'%d \n',[1,1]); [~] = fgetl(fid);
prm.rwd = fscanf(fid,'%f \n',[1,1]); [~] = fgetl(fid);
prm.nrep= fscanf(fid,'%i \n',[1,1]); [~] = fgetl(fid);
prm.efrq= fscanf(fid,'%i \n',[1,1]); [~] = fgetl(fid);
prm.tmax= fscanf(fid,'%f \n',[1,1]);
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
fprintf('File Rigid boundary       : %s \n',prm.filebound)
fprintf('File Internal deformation : %s \n',prm.fileinternal)
fprintf('File UDline settings      : %s \n',prm.fileUDconfig)
fprintf('DIRResult                 : %s \n',prm.dirresult)
fprintf('GPUdev (CPU:99)           : %i \n',prm.gpu)
fprintf('ITR(Max_Nitr)             : %i \n',prm.itr)
fprintf('ITR(Threshold_Nitr)       : %i \n',prm.thr)
fprintf('CHA(Chain)                : %i \n',prm.cha)
fprintf('KEP(KEEP)                 : %i \n',prm.kep)
fprintf('RWD(Walk_dis)             : %4.2f \n',prm.rwd)
fprintf('Number of Replica         : %i \n',prm.nrep)
fprintf('Exchange Frequency        : %i \n',prm.efrq)
fprintf('Maximum temperature       : %5.1f \n',prm.tmax) 
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
    if (lca(1)==true && lca(end)==true) && (lca(2)==false && lca(end-1)==false)  % eliminate Quadruple-junction
      lca( 1 ) = false;
      lca(end) = false;
    end
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
function [blk] = ReadBlockInterface(blk,prm,obs)
% Coded by Takeo Ito 2016/12/21 (ver 1.0)
%
int_tri       =   50;
dep_limit     = -100;
dep_limit_low =  -20;
dirblk        = prm.dirblock_interface;
blk(1).f     = [];
blk(1).da    = [];
blk(1).nv    = [];
blk(1).st    = [];
blk(1).dp    = [];
blk(1).phi   = [];
blk(1).theta = [];
blk(1).clon  = [];
blk(1).clat  = [];
blk(1).cdep  = [];
blk(1).nt    =  0; % Number of trimeshes
blk(1).ntmec =  0;
blk(1).ntkin =  0;
alat = mean(obs(1).alat(:));
alon = mean(obs(1).alon(:));
% 
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    blk(1).bound(nb1,nb2).flag1 = 0;
    pre_tri_f = fullfile(dirblk,['triB_',num2str(nb1),'_',num2str(nb2),'.txt']); 
    fid = fopen(pre_tri_f,'r');
    if fid >= 0
      fprintf('block interface: %2i  %2i \n',nb1,nb2)
      fprintf('read interface tri boudary file : %s \n',pre_tri_f)
      for ii = 1:size(blk(1).dipbo,1)
        dippingid = ismember([nb1 nb2],blk(1).dipbo(ii,2:3));
        ispair    = sum(dippingid);
        if ispair == 2 && ismember(blk(1).dipbo(ii,1),[0:2])
          blk(1).bound(nb1,nb2).flag1 = blk(1).dipbo(ii,1);
          break
        end
      end
      nf   = 0;
      blon = zeros(1,3);
      blat = zeros(1,3);
      bdep = zeros(1,3);
      type = 0;
      while 1
        nf    = nf+1;
        tline = fgetl(fid); if ~ischar(tline); break; end
        lchar = strsplit(tline); if strcmpi(lchar{1},''); break; end
        loc_f = fscanf(fid,'%f %f %f \n', [3 3]);
        blon(nf,:) = loc_f(1,:);  % lon
        blat(nf,:) = loc_f(2,:);  % lat
        bdep(nf,:) = loc_f(3,:);  % hight
        if size(lchar,2) < 2
          type(nf) = blk(1).bound(nb1,nb2).flag1;
        else
          type(nf) = str2num(lchar{2});
        end
        tline = fgetl(fid); if ~ischar(tline); break; end
        lchar = strsplit(tline); if strcmpi(lchar{1},''); break; end
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
        type = type(id);
      end
      blk(1).bound(nb1,nb2).blon = blon;  % lon
      blk(1).bound(nb1,nb2).blat = blat;  % lat
      blk(1).bound(nb1,nb2).bdep = bdep;  % hight
      blk(1).bound(nb1,nb2).type = type;  % mesh type
    else
      sub_f = fullfile(dirblk,['B_',num2str(nb1),'_',num2str(nb2),'.txt']);
      fid   = fopen(sub_f,'r');
      if fid >= 0
        fprintf('block interface: %2i  %2i \n',nb1,nb2)
        fprintf('read interface boudary shape file : %s \n',sub_f)
        for ii = 1:size(blk(1).dipbo,1)
          dippingid = ismember([nb1 nb2],blk(1).dipbo(ii,2:3));
          ispair    = sum(dippingid);
          if ispair == 2 && ismember(blk(1).dipbo(ii,1),[0:2])
            blk(1).bound(nb1,nb2).flag1 = blk(1).dipbo(ii,1);
            break
          end
        end
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
        blk(1).bound(nb1,nb2).type = blk(1).bound(nb1,nb2).flag1.*ones(1,size(bstri,1));
      end
    end

    if size(blk(1).bound(nb1,nb2).blon,1) ~= 0
      % Save trimeshes to text file
      out_tri_f = fullfile(prm.dirblock_interface,['triB_',num2str(nb1),'_',num2str(nb2),'.out']);
      nlen      = length(blk(1).bound(nb1,nb2).blat(:,1));
      fid_out   = fopen(out_tri_f,'wt');
      for ntri = 1:nlen
        fprintf(fid_out,'> %i\n',blk(1).bound(nb1,nb2).type(ntri));
        fprintf(fid_out,'%f %f %f \n%f %f %f \n%f %f %f \n%f %f %f \n',...
                [blk(1).bound(nb1,nb2).blon(ntri,1), blk(1).bound(nb1,nb2).blat(ntri,1), blk(1).bound(nb1,nb2).bdep(ntri,1);...
                 blk(1).bound(nb1,nb2).blon(ntri,2), blk(1).bound(nb1,nb2).blat(ntri,2), blk(1).bound(nb1,nb2).bdep(ntri,2);...
                 blk(1).bound(nb1,nb2).blon(ntri,3), blk(1).bound(nb1,nb2).blat(ntri,3), blk(1).bound(nb1,nb2).bdep(ntri,3);...
                 blk(1).bound(nb1,nb2).blon(ntri,1), blk(1).bound(nb1,nb2).blat(ntri,1), blk(1).bound(nb1,nb2).bdep(ntri,1)]');
      end
      fclose(fid_out);

      % Calc angle and area of trimeshes
      for nf=1:size(blk(1).bound(nb1,nb2).blon,1)
        [trix, triy] = PLTXY(blk(1).bound(nb1,nb2).blat(nf,:), blk(1).bound(nb1,nb2).blon(nf,:), alat, alon);
        triz         = -1.*blk(1).bound(nb1,nb2).bdep(nf,:);
        f_loc       = [trix; triy; triz];
        [f,da,nv,st,dp,phi,theta] = EstFaultTri(f_loc);
        blk(1).bound(nb1,nb2).f(    nf,:) =     f ;
        blk(1).bound(nb1,nb2).da(   nf,:) =    da ;
        blk(1).bound(nb1,nb2).nv(   nf,:) =    nv';
        blk(1).bound(nb1,nb2).st(   nf,:) =    st ;
        blk(1).bound(nb1,nb2).dp(   nf,:) =    dp ;
        blk(1).bound(nb1,nb2).phi(  nf,:) =   phi ;
        blk(1).bound(nb1,nb2).theta(nf,:) = theta ;
      end
      
      % Stack calculated value
      blk(1).f     = [blk(1).f    ; blk(1).bound(nb1,nb2).f           ];
      blk(1).da    = [blk(1).da   ; blk(1).bound(nb1,nb2).da          ];
      blk(1).nv    = [blk(1).nv   ; blk(1).bound(nb1,nb2).nv          ];
      blk(1).st    = [blk(1).st   ; blk(1).bound(nb1,nb2).st          ];
      blk(1).dp    = [blk(1).dp   ; blk(1).bound(nb1,nb2).dp          ];
      blk(1).phi   = [blk(1).phi  ; blk(1).bound(nb1,nb2).phi         ];
      blk(1).theta = [blk(1).theta; blk(1).bound(nb1,nb2).theta       ];
      blk(1).clon  = [blk(1).clon ; mean(blk(1).bound(nb1,nb2).blon,2)];
      blk(1).clat  = [blk(1).clat ; mean(blk(1).bound(nb1,nb2).blat,2)];
      blk(1).cdep  = [blk(1).cdep ; mean(blk(1).bound(nb1,nb2).bdep,2)];

      % Grouping boundary type
      blk(1).bound(nb1,nb2).flag1 = 0;
      for ii = 1:size(blk(1).dipbo,1)
        dippingid = ismember([nb1 nb2],blk(1).dipbo(ii,2:3));
        ispair    = sum(dippingid);
        if ispair == 2 && ismember(blk(1).dipbo(ii,1),[0:2])
          blk(1).bound(nb1,nb2).flag1 = blk(1).dipbo(ii,1);
          break
        end
      end
      blk(1).bound(nb1,nb2).flag2 = 0;
      for ii = 1:size(blk(1).mechbo,1)
        mechcpid  = ismember([nb1 nb2],blk(1).mechbo(ii,:));
        ispair    = sum(mechcpid);
        if ispair == 2
          blk(1).bound(nb1,nb2).flag2 = 1;
          break
        end
      end

      blk(1).nt = blk(1).nt + size(blk(1).bound(nb1,nb2).blon,1);
      if blk(1).bound(nb1,nb2).flag2 == 1
        blk(1).ntmec = blk(1).ntmec + size(blk(1).bound(nb1,nb2).blon,1);
      else
        blk(1).ntkin = blk(1).ntkin + size(blk(1).bound(nb1,nb2).blon,1);
      end
    end
  end
end
end

%% Read locked patches
function [blk] = ReadLockedPatch(blk,prm)
%    Test version coded by H. Kimura 2019/1/29
% Revised version coded by H. Kimura 2019/2/5

for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    if blk(1).bound(nb1,nb2).flag2 == 1
      blk(1).bound(nb1,nb2).patchid = 0;
      patchfile = fullfile(prm.dirblock_interface,['patchb_',num2str(nb1),'_',num2str(nb2),'.txt']);
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
      else
        error(['Not found', patchfile]);
      end
    end
  end
end

fprintf('=== Read Locked Patches=== \n');
end

%% Show block map
function ShowBlockBound(blk)
%
minlat  = []; minlon  = []; maxlat  = []; maxlon  = [];
minlatc = []; minlonc = []; maxlatc = []; maxlonc = [];
for nb = 1:blk(1).nblock
  minlatc = min([(min(blk(nb).lat) + max(blk(nb).lat))/2;minlatc]); 
  minlonc = min([(min(blk(nb).lon) + max(blk(nb).lon))/2;minlonc]); 
  maxlatc = max([(min(blk(nb).lat) + max(blk(nb).lat))/2;maxlatc]); 
  maxlonc = max([(min(blk(nb).lon) + max(blk(nb).lon))/2;maxlonc]); 
  minlat = min([blk(nb).lat;minlat]); 
  minlon = min([blk(nb).lon;minlon]); 
  maxlat = max([blk(nb).lat;maxlat]); 
  maxlon = max([blk(nb).lon;maxlon]); 
end

f300 = figure(300); clf(300)
f300.Position = [60,100,1400,800];
% Blocks (world wide)
figure(300); subplot(1,2,1)
h = worldmap([minlat,maxlat],[minlon,maxlon]);
getm(h, 'MapProjection');
geoshow('landareas.shp', 'FaceColor', [0.15 0.5 0.15])
for nb = 1:blk(1).nblock
  hold on
  plotm(blk(nb).lat,blk(nb).lon,'red')
  hold on
  textm(mean(blk(nb).lat),mean(blk(nb).lon),int2str(nb))
end
drawnow

figure(300); subplot(1,2,2)
% Blocks (regional)
h = worldmap([minlatc,maxlatc],[minlonc,maxlonc]);
getm(h, 'MapProjection');
geoshow('landareas.shp', 'FaceColor', [0.15 0.5 0.15])
for nb1=1:blk(1).nblock
  hold on
  plotm(blk(nb1).lat,blk(nb1).lon,'red')
  hold on
  textm(mean(blk(nb).lat),mean(blk(nb).lon),int2str(nb))
  for nb2=nb1+1:blk(1).nblock
    if ~isempty(blk(1).bound(nb1,nb2).lat)
      hold on
      plotm(blk(1).bound(nb1,nb2).lat,blk(1).bound(nb1,nb2).lon,'o')
    end
  end
end
drawnow

% Block and trimesh
figure(310); clf(310)
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    nf = size(blk(1).bound(nb1,nb2).blon,1);
    if nf ~= 0
      patch(blk(1).bound(nb1,nb2).blon',blk(1).bound(nb1,nb2).blat',blk(1).bound(nb1,nb2).bdep',blk(1).bound(nb1,nb2).bdep');
      hold on
    end
  end
end
drawnow 
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
savefig(300,fullfile(f_fir,'region_map'))
savefig(310,fullfile(f_fir,'mesh_map'))
if ~strcmpi(prm.method,'Forward')
  movefile([prm.dirresult,'/cha_test*.mat'],a_dir)
  savefig(100,fullfile(f_fir,'coupling_bslip'))
  savefig(130,fullfile(f_fir,'pole'))
  savefig(140,fullfile(f_fir,'vector'))
  savefig(160,fullfile(f_fir,'strain'))
else
  savefig(100,fullfile(f_fir,'bslip_vector'))
end
% 
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
  fprintf(       'BLOCK=%2d NUM_OBS=%2d Sigma^2=%5.2f ',n,obs(n).nblk,sig)
  fprintf(logfid,'BLOCK=%2d NUM_OBS=%2d Sigma^2=%5.2f ',n,obs(n).nblk,sig);
  [latp,lonp,ang] = xyzp2lla(pole(1),pole(2),pole(3));
  fprintf(       'Lat:%7.2f deg. Lon:%8.2f deg. Ang:%9.2e deg./m.y. \n',latp,lonp,ang);    
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

%% Read udline share and interp interval
function [blk] = ReadUDlineSettings(blk,prm)
%--example from here--
% # Interp intervals
% 5 8 2
% 6 8 2
% 7 8 4
% 10 10 3
% 11 13 5
% # UDline share
% 5 8  6 8  7 8
% 10 13  11 13
% --- END HERE ---
%--end of example--
fid = fopen(prm.fileUDconfig,'r');
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    blk(1).bound(nb1,nb2).udshare  = [];
    blk(1).bound(nb1,nb2).udinterp = [];
  end
end
if fid ~= 0
  tline = char(fgetl(fid));
  while 1
    switch tline
      case '# UDline share'
        while 1
          tline = char(fgetl(fid));
          Tline = strtrim(strsplit(tline));
          if ~or(strcmpi(Tline(1),'---'),or(strcmpi(Tline(1),'#'),strcmpi(Tline(1),'')))
            Tline=str2num(char(Tline))';
            if ~isempty(Tline)
              for n = 2:size(Tline,2)/2
                blk(1).bound(Tline(2*n-1),Tline(2*n)).udshare = Tline(1:2);
              end
            end
          else
            break
          end
        end
      case '# Interp intervals'
        while 1
          tline = char(fgetl(fid));
          Tline = strtrim(strsplit(tline));
          if ~or(strcmpi(Tline(1),'---'),or(strcmpi(Tline(1),'#'),strcmpi(Tline(1),'')))
            Tline=str2num(char(Tline));
            if ~isempty(Tline)
              blk(1).bound(Tline(1),Tline(2)).udinterp = Tline(3);
            end
          else
            break
          end
        end
      otherwise
        tline = char(fgetl(fid));
    end
    if strcmpi(tline,'--- END HERE ---'); break; end
  end
end

end

%% Read dipping block boundary
function [blk,prm] = ReadBoundaryTypes(blk,prm)
% Note:
% Prepare the export parameter file in the 'PARAMETER' folder as bellows,
%--example from here--
% # Dipping boundaries
% 8 7 8
% 11 8 11
% 9 8 9
% 10 9 10
% 12 10 12
% # Mechanical coupled boundaries
% 8 11
% 8 9
% 9 10
% 10 12
% --- END HERE ---
%--end of example--

fid = fopen(prm.filebound,'r');
blk(1).dipbo  = [];
blk(1).mechbo = [];
if fid ~= 0
  tline = char(fgetl(fid));
  while 1
    switch tline
      case '# Dipping boundaries'
        while 1
          tline = char(fgetl(fid));
          Tline = strtrim(strsplit(tline));
          if ~or(strcmpi(Tline(1),'---'),or(strcmpi(Tline(1),'#'),strcmpi(Tline(1),'')))
            Tline=str2num(char(Tline));
            if ~isempty(Tline)
              blk(1).dipbo = [blk(1).dipbo; Tline'];
            end
          else
            break
          end
        end
      case '# Mechanical coupled boundaries'
        while 1
          tline = char(fgetl(fid));
          Tline = strtrim(strsplit(tline));
          if ~or(strcmpi(Tline(1),'---'),or(strcmpi(Tline(1),'#'),strcmpi(Tline(1),'')))
            Tline=str2num(char(Tline));
            if ~isempty(Tline)
              blk(1).mechbo = [blk(1).mechbo; Tline'];
            end
          else
            break
          end
        end
      otherwise
        tline = char(fgetl(fid));
    end
    if strcmpi(tline,'--- END HERE ---'); break; end
  end
end

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
function [blk,tri] = GreenTri(blk,obs,prm)
% Coded by Hiroshi Kimura 2019/01/28 (ver 1.0)
pr = 0.25;
nd = size(obs(1).alat,2);
nfall = blk(1).nt;
c1=cos(blk(1).phi);
s1=sin(blk(1).phi);
c2=cos(-blk(1).theta);
s2=sin(-blk(1).theta);
%
alat = mean(obs(1).alat(:));
alon = mean(obs(1).alon(:));
[obsx,obsy] = PLTXY(obs(1).alat,obs(1).alon,alat,alon);
obsz = -1e-3.*obs(1).ahig;
%
[tricx,tricy] = PLTXY(blk(1).clat,blk(1).clon,alat,alon);
tricz = -blk(1).cdep;
% 
mc           = 1;
tri(1).obsdis= [];
blk(1).nb    = 0;
blk(1).nbmec = 0;
blk(1).nbkin = 0;
tri(1).idstr = false(3*nfall,1);
tri(1).iddip = false(3*nfall,1);
tri(1).idtns = false(3*nfall,1);
for nb1 = 1:blk(1).nblock
  [tri(nb1).localx, tri(nb1).localy] = PLTXY(blk(nb1).lat,blk(nb1).lon,alat,alon);
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
        tri(1).bound(nb1,nb2).gsstr = zeros(6*nfall,nf);
        tri(1).bound(nb1,nb2).gsdip = zeros(6*nfall,nf);
        tri(1).bound(nb1,nb2).gstns = zeros(6*nfall,nf);
        tri(1).bound(nb1,nb2).cf   =  ones(3*nf,1);
        tri(1).bound(nb1,nb2).inv  = zeros(3*nf,1);
        %
        fprintf('==================\n Block %2i : Block %2i \n Number of TRI sub-faults : %4i \n',nb1,nb2,nf)
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
          % Green's function relates to strike slip
          u = CalcTriDisps(obsx,obsy,obsz,trix,triy,triz,pr,1,0,0);
          s = CalcTriStrains(tricx,tricy,tricz,trix,triy,triz,pr,1,0,0);
          tri(1).bound(nb1,nb2).gustr(1:3:3*nd   ,n) =  u.x;   % E
          tri(1).bound(nb1,nb2).gustr(2:3:3*nd   ,n) =  u.y;   % N
          tri(1).bound(nb1,nb2).gustr(3:3:3*nd   ,n) = -u.z;   % U
          tri(1).bound(nb1,nb2).gsstr(1:6:6*nfall,n) =  s.xx;  % exx
          tri(1).bound(nb1,nb2).gsstr(2:6:6*nfall,n) =  s.xy;  % exy
          tri(1).bound(nb1,nb2).gsstr(3:6:6*nfall,n) = -s.xz;  % exz
          tri(1).bound(nb1,nb2).gsstr(4:6:6*nfall,n) =  s.yy;  % eyy
          tri(1).bound(nb1,nb2).gsstr(5:6:6*nfall,n) = -s.yz;  % eyz
          tri(1).bound(nb1,nb2).gsstr(6:6:6*nfall,n) = -s.zz;  % ezz
          % Green's function relates to tensile slip
          u = CalcTriDisps(obsx,obsy,obsz,trix,triy,triz,pr,0,1,0);
          s = CalcTriStrains(tricx,tricy,tricz,trix,triy,triz,pr,0,1,0);
          tri(1).bound(nb1,nb2).gutns(1:3:3*nd   ,n) =  u.x;   % E
          tri(1).bound(nb1,nb2).gutns(2:3:3*nd   ,n) =  u.y;   % N
          tri(1).bound(nb1,nb2).gutns(3:3:3*nd   ,n) = -u.z;   % U 
          tri(1).bound(nb1,nb2).gstns(1:6:6*nfall,n) =  s.xx;  % exx
          tri(1).bound(nb1,nb2).gstns(2:6:6*nfall,n) =  s.xy;  % exy
          tri(1).bound(nb1,nb2).gstns(3:6:6*nfall,n) = -s.xz;  % exz
          tri(1).bound(nb1,nb2).gstns(4:6:6*nfall,n) =  s.yy;  % eyy
          tri(1).bound(nb1,nb2).gstns(5:6:6*nfall,n) = -s.yz;  % eyz
          tri(1).bound(nb1,nb2).gstns(6:6:6*nfall,n) = -s.zz;  % ezz
          % Green's function relates to dip slip
          u = CalcTriDisps(obsx,obsy,obsz,trix,triy,triz,pr,0,0,1);
          s = CalcTriStrains(tricx,tricy,tricz,trix,triy,triz,pr,0,0,1);
          tri(1).bound(nb1,nb2).gudip(1:3:3*nd   ,n) =  u.x;   % E
          tri(1).bound(nb1,nb2).gudip(2:3:3*nd   ,n) =  u.y;   % N
          tri(1).bound(nb1,nb2).gudip(3:3:3*nd   ,n) = -u.z;   % U
          tri(1).bound(nb1,nb2).gsdip(1:6:6*nfall,n) =  s.xx;  % exx
          tri(1).bound(nb1,nb2).gsdip(2:6:6*nfall,n) =  s.xy;  % exy
          tri(1).bound(nb1,nb2).gsdip(3:6:6*nfall,n) = -s.xz;  % exz
          tri(1).bound(nb1,nb2).gsdip(4:6:6*nfall,n) =  s.yy;  % eyy
          tri(1).bound(nb1,nb2).gsdip(5:6:6*nfall,n) = -s.yz;  % eyz
          tri(1).bound(nb1,nb2).gsdip(6:6:6*nfall,n) = -s.zz;  % ezz
          if mod(n,ceil(nf/3)) == 1
            fprintf('MAKE GREEN at TRI sub-faults : %4i / %4i \n',n,nf)
          end
          [tri] = CorrectFactor(blk,tri,nb1,nb2,dp,n,nf);
          [tri] = DiscriminateDirection(blk,tri,nb1,nb2,trix,triy,n,nf);
          [tri] = trans_xyz2strdip(blk,tri,nb1,nb2,n,nfall,c1,s1,c2,s2);
        end
        trisave = tri(1).bound(nb1,nb2);
        save(trimat,'trisave','-v7.3');
      end
      tri(1).idstr(mc     :mc+  nf-1) = true;
      tri(1).iddip(mc+  nf:mc+2*nf-1) = true;
      tri(1).idtns(mc+2*nf:mc+3*nf-1) = true;
      blk(1).nb = blk(1).nb + 1;
      if blk(1).bound(nb1,nb2).flag2 == 1
        blk(1).nbmec = blk(1).nbmec + 1;
      else
        blk(1).nbkin = blk(1).nbkin + 1;
      end
      mc = mc + 3*nf;
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
switch blk(1).bound(nb1,nb2).type(nf)
  case {1,2}
    cf = 1/sqrt(dp(1)^2+dp(2)^2);      % 1=sqrt(dp(1)^2+dp(2)^2+dp(3)^2): norm of dp
    if cf==Inf, cf=1; end
    tri(1).bound(nb1,nb2).cf(nf+n) = cf;
  case 0
    return
end
end

%% DISCRIMINATE BOUNDARY TYPE AND SUBFAULT SURFACE DIRECTION
function [tri] = DiscriminateDirection(blk,tri,nb1,nb2,trix,triy,n,nf)
% Coded by H.Kimura 2017/4/28 (test ver.)
% Modified by H.Kimura 2018/2/6
switch blk(1).bound(nb1,nb2).type(nf)
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
    [in,on] = inpolygon(trixc,triyc,tri(nb1).localx,tri(nb1).localy);
    if in==1 && on~=1
      tri(1).bound(nb1,nb2).inv(     n) = 1;
      tri(1).bound(nb1,nb2).inv(  nf+n) = 0;
      tri(1).bound(nb1,nb2).inv(2*nf+n) = 1;
    elseif in==1 && on==1
      tric = [trixc triyc 0];
      uv   = [0 0 1];
      nv   = cross(uv,tri(1).bound(nb1,nb2).st(n,:));
      cnv  = tric+nv;
      if inpolygon(cnv(1),cnv(2),tri(nb1).localx,tri(nb1).localy)==1
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
function [tri] = trans_xyz2strdip(blk,tri,nb1,nb2,n,nfall,c1,s1,c2,s2)
% This function transforms a strain tensor from xyz to fault strike-dip.
% Strain of strike direction on the fault corresponds to ezx',
% strain of dip direction on the fault corresponds to ezy', and
% strain of tensile direction on the fault corresponds to ezz'.

% c1=cos(tri(1).bound(nb1,nb2).phi(n));
% s1=sin(tri(1).bound(nb1,nb2).phi(n));
% c2=cos(-tri(1).bound(nb1,nb2).theta(n));
% s2=sin(-tri(1).bound(nb1,nb2).theta(n));

[sst,sdp,sts] = calctrans(tri(1).bound(nb1,nb2).gsstr(:,n),c1,s1,c2,s2);  % response to strike slip
tri(1).bound(nb1,nb2).gsstrT(1:3:3*nfall,n)=sst;  % str
tri(1).bound(nb1,nb2).gsstrT(2:3:3*nfall,n)=sdp;  % dip
tri(1).bound(nb1,nb2).gsstrT(3:3:3*nfall,n)=sts;  % tns
[sst,sdp,sts] = calctrans(tri(1).bound(nb1,nb2).gstns(:,n),c1,s1,c2,s2);  % response to tensile slip
tri(1).bound(nb1,nb2).gstnsT(1:3:3*nfall,n)=sst;  % str
tri(1).bound(nb1,nb2).gstnsT(2:3:3*nfall,n)=sdp;  % dip
tri(1).bound(nb1,nb2).gstnsT(3:3:3*nfall,n)=sts;  % tns
[sst,sdp,sts] = calctrans(tri(1).bound(nb1,nb2).gsdip(:,n),c1,s1,c2,s2);  % response to dip slip
tri(1).bound(nb1,nb2).gsdipT(1:3:3*nfall,n)=sst;  % str
tri(1).bound(nb1,nb2).gsdipT(2:3:3*nfall,n)=sdp;  % dip
tri(1).bound(nb1,nb2).gsdipT(3:3:3*nfall,n)=sts;  % tns

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
d(1).ind   = find(tmp.err ~= 0)';
d(1).obs   = tmp.obs(d(1).ind)';
d(1).err   = tmp.err(d(1).ind)';
d(1).mcid  = zeros(3*blk(1).ntkin,blk(1).nbkin);  % kinematic couple
d(1).maid  = zeros(3*blk(1).ntmec,blk(1).nbmec);  % mechanical couple
d(1).idmec = false(3*blk(1).nt   ,           1);
%
% (G(1).C * (( G(1).T * ( G(1).B1 - G(1).B2 ) * Mp)*Mc ) + G(1).P * Mp
% 
% (G(1).C * (( G(1).T * ( G(1).B1 - G(1).B2 ) * Mp)*Mc ) + G(1).P * Mp +
% G(1).I * Mi  % including internal deformation
%
G(1).t  = zeros(3*blk(1).nt,2.*blk(1).nt);
G(1).b  = zeros(2*blk(1).nt,3.*blk(1).nblock);
G(1).zc = zeros(blk(1).ntmec,1);
tmp.p   = zeros(3*nobs,3.*blk(1).nblock);
tmp.i   = zeros(3*nobs,3.*blk(1).nblock);
tmp.s   = zeros(3*blk(1).nt,3*blk(1).nt);
tmp.cf  = ones(3*blk(1).nt,1);
tmp.inv = zeros(3*blk(1).nt,1);
tmp.zlimit = zeros(blk(1).naspline,2);
% 
nasp = blk(1).naspline;
% 
mf1 = 1;
mf2 = 1;
mf3 = 1;
ma  = 1;
ms  = 1;
mm1 = 1;
mm3 = 1;
mk1 = 1;
mk3 = 1;
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    nf = size(tri(1).bound(nb1,nb2).clon,2);
    if nf ~= 0
      if blk(1).bound(nb1,nb2).flag2 == 1
        % d(1).maid(mm3:mm3+3*nf-1,mm1:mm1+nf-1) = repmat(eye(nf),3,1);
        d(1).maid(mm3:mm3+3*nf-1,mm1:mm1+nf-1) = [repmat(eye(nf),2,1); zeros(nf)];  % due to "zero tns-strain"
        nint = blk(1).bound(nb1,nb2).interpint;
        nasp = blk(1).bound(nb1,nb2).naspline ;
        nstp = nint * (nasp - 1) + 1;
        d(1).idmec(mf3:mf3+3*nf-1) = true;
        G(1).zc(mm1:mm1+nf-1) = -blk(1).cdep(mf1:mf1+nf-1);
        if size(blk(1).bound(nb1,nb2).udshare,1) ~= 0
          matmp = tmp.pair(blk(1).bound(nb1,nb2).udshare(1),blk(1).bound(nb1,nb2).udshare(2)).matmp;
          tmp.idMit(ms:ms+nstp-1,matmp:matmp+nasp-1) = blk(1).bound(nb1,nb2).interpid;
        else
          tmp.zlimit(ma:ma+nasp-1,1) = blk(1).bound(nb1,nb2).asp_depd;
          tmp.zlimit(ma:ma+nasp-1,2) = blk(1).bound(nb1,nb2).asp_depu;
          tmp.idMit(ms:ms+nstp-1,   ma:   ma+nasp-1) = blk(1).bound(nb1,nb2).interpid;
          tmp.pair(nb1,nb2).matmp = ma;
          ma = ma + nasp;
        end
        tmp.idstrip(mm1:mm1+nf-1,ms:ms+nstp-1) = blk(1).bound(nb1,nb2).stripid ;
        mm1 = mm1 +   nf;
        mm3 = mm3 + 3*nf;
        ms = ms + nstp;
      else
        d(1).mcid(mk3:mk3+3*nf-1,mk1:mk1+nf-1) = repmat(eye(nf),3,1);
        mk1 = mk1 +   nf;
        mk3 = mk3 + 3*nf;
      end
      tmp.c(1:3*nobs,mf3     :mf3+  nf-1) = tri(1).bound(nb1,nb2).gustr;
      tmp.c(1:3*nobs,mf3+  nf:mf3+2*nf-1) = tri(1).bound(nb1,nb2).gudip;
      tmp.c(1:3*nobs,mf3+2*nf:mf3+3*nf-1) = tri(1).bound(nb1,nb2).gutns;
      tmp.cf( mf3   :mf3+3*nf-1) = tri(1).bound(nb1,nb2).cf ;
      tmp.inv(mf3   :mf3+3*nf-1) = tri(1).bound(nb1,nb2).inv;
      G(1).t(mf3   :mf3+  nf-1,mf2   :mf2+  nf-1) = diag(tri(1).bound(nb1,nb2).st(:,1));
      G(1).t(mf3+nf:mf3+2*nf-1,mf2   :mf2+  nf-1) = diag(tri(1).bound(nb1,nb2).dp(:,1));
      G(1).t(mf3   :mf3+  nf-1,mf2+nf:mf2+2*nf-1) = diag(tri(1).bound(nb1,nb2).st(:,2));
      G(1).t(mf3+nf:mf3+2*nf-1,mf2+nf:mf2+2*nf-1) = diag(tri(1).bound(nb1,nb2).dp(:,2));
      G(1).b(mf2   :mf2+  nf-1,3*nb1-2) = -1.*(-tri(1).bound(nb1,nb2).oxyz(:,7).*tri(1).bound(nb1,nb2).oxyz(:,3));
      G(1).b(mf2   :mf2+  nf-1,3*nb1-1) = -1.*(-tri(1).bound(nb1,nb2).oxyz(:,5).*tri(1).bound(nb1,nb2).oxyz(:,3));
      G(1).b(mf2   :mf2+  nf-1,3*nb1  ) = -1.*( tri(1).bound(nb1,nb2).oxyz(:,5).*tri(1).bound(nb1,nb2).oxyz(:,2)...
                                      +tri(1).bound(nb1,nb2).oxyz(:,7).*tri(1).bound(nb1,nb2).oxyz(:,1));
      G(1).b(mf2+nf:mf2+2*nf-1,3*nb1-2) = -1.*( tri(1).bound(nb1,nb2).oxyz(:,4).*tri(1).bound(nb1,nb2).oxyz(:,5).*tri(1).bound(nb1,nb2).oxyz(:,3)...
                                      +tri(1).bound(nb1,nb2).oxyz(:,6).*tri(1).bound(nb1,nb2).oxyz(:,2));
      G(1).b(mf2+nf:mf2+2*nf-1,3*nb1-1) = -1.*(-tri(1).bound(nb1,nb2).oxyz(:,4).*tri(1).bound(nb1,nb2).oxyz(:,7).*tri(1).bound(nb1,nb2).oxyz(:,3)...
                                      -tri(1).bound(nb1,nb2).oxyz(:,6).*tri(1).bound(nb1,nb2).oxyz(:,1));
      G(1).b(mf2+nf:mf2+2*nf-1,3*nb1  ) = -1.*( tri(1).bound(nb1,nb2).oxyz(:,4).*tri(1).bound(nb1,nb2).oxyz(:,7).*tri(1).bound(nb1,nb2).oxyz(:,2)...
                                      -tri(1).bound(nb1,nb2).oxyz(:,4).*tri(1).bound(nb1,nb2).oxyz(:,5).*tri(1).bound(nb1,nb2).oxyz(:,1));
      G(1).b(mf2   :mf2+  nf-1,3*nb2-2) = -G(1).b(mf2   :mf2+  nf-1,3*nb1-2);
      G(1).b(mf2   :mf2+  nf-1,3*nb2-1) = -G(1).b(mf2   :mf2+  nf-1,3*nb1-1);
      G(1).b(mf2   :mf2+  nf-1,3*nb2  ) = -G(1).b(mf2   :mf2+  nf-1,3*nb1  );
      G(1).b(mf2+nf:mf2+2*nf-1,3*nb2-2) = -G(1).b(mf2+nf:mf2+2*nf-1,3*nb1-2);
      G(1).b(mf2+nf:mf2+2*nf-1,3*nb2-1) = -G(1).b(mf2+nf:mf2+2*nf-1,3*nb1-1);
      G(1).b(mf2+nf:mf2+2*nf-1,3*nb2  ) = -G(1).b(mf2+nf:mf2+2*nf-1,3*nb1  );
      tmp.s(tri(1).idstr,mf3     :mf3+  nf-1) = tri(1).bound(nb1,nb2).gsstrT(1:3:end,:);
      tmp.s(tri(1).idstr,mf3+  nf:mf3+2*nf-1) = tri(1).bound(nb1,nb2).gsdipT(1:3:end,:);
      %       tmp.s(tri(1).idstr,mc+2*nf:mc+3*nf-1) = tri(1).bound(nb1,nb2).gstnsT(1:3:end,:);
      tmp.s(tri(1).idstr,mf3+2*nf:mf3+3*nf-1) = zeros(size(tri(1).bound(nb1,nb2).gstnsT(1:3:end,:)));
      tmp.s(tri(1).iddip,mf3     :mf3+  nf-1) = tri(1).bound(nb1,nb2).gsstrT(2:3:end,:);
      tmp.s(tri(1).iddip,mf3+  nf:mf3+2*nf-1) = tri(1).bound(nb1,nb2).gsdipT(2:3:end,:);
      %       tmp.s(tri(1).iddip,mc+2*nf:mc+3*nf-1) = tri(1).bound(nb1,nb2).gstnsT(2:3:end,:);
      tmp.s(tri(1).iddip,mf3+2*nf:mf3+3*nf-1) = zeros(size(tri(1).bound(nb1,nb2).gstnsT(2:3:end,:)));
      %       tmp.s(tri(1).idtns,mc     :mc+  nf-1) = tri(1).bound(nb1,nb2).gsstrT(3:3:end,:);
      %       tmp.s(tri(1).idtns,mc+  nf:mc+2*nf-1) = tri(1).bound(nb1,nb2).gsdipT(3:3:end,:);
      %       tmp.s(tri(1).idtns,mc+2*nf:mc+3*nf-1) = tri(1).bound(nb1,nb2).gstnsT(3:3:end,:);
      tmp.s(tri(1).idtns,mf3     :mf3+  nf-1) = zeros(size(tri(1).bound(nb1,nb2).gsstrT(3:3:end,:)));
      tmp.s(tri(1).idtns,mf3+  nf:mf3+2*nf-1) = zeros(size(tri(1).bound(nb1,nb2).gsdipT(3:3:end,:)));
      tmp.s(tri(1).idtns,mf3+2*nf:mf3+3*nf-1) = zeros(size(tri(1).bound(nb1,nb2).gstnsT(3:3:end,:)));
      %
      mf1 = mf1 +   nf;
      mf2 = mf2 + 2*nf;
      mf3 = mf3 + 3*nf;
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
tmp.c     = tmp.c(d(1).ind,:);
tmp.tb    = G(1).t*G(1).b;
tmp.cfinv = tmp.cf.*tmp.inv;
tmp.zu     = [zeros(blk(1).naspline),   eye(blk(1).naspline)];
tmp.zd     = [  eye(blk(1).naspline), zeros(blk(1).naspline)];
tmp.zlimit = reshape(tmp.zlimit, 2*blk(1).naspline, 1);
G(1).tb_kin = sparse(tmp.tb(~d(1).idmec,:));
G(1).tb_mec = sparse(tmp.tb( d(1).idmec,:));
G(1).c_kin  = tmp.c(:,~d(1).idmec);
G(1).c_mec  = tmp.c(:, d(1).idmec);
G(1).p      = tmp.p(d(1).ind,:);
G(1).i      = tmp.i(d(1).ind,:);
% G(1).s      = tmp.s(d(1).idmec,d(1).idmec);
tmpGs = tmp.s(d(1).idmec,d(1).idmec);
tmpGs(abs(tmpGs) < max(abs(tmpGs)).*(1e-2/blk(1).ntmec)) = 0;
G(1).s      = tmpGs;
G(1).zu     = tmp.idstrip * tmp.idMit * tmp.zu;
G(1).zd     = tmp.idstrip * tmp.idMit * tmp.zd;
G(1).zulim  = tmp.idstrip * tmp.idMit * tmp.zu * tmp.zlimit;
G(1).zdlim  = tmp.idstrip * tmp.idMit * tmp.zd * tmp.zlimit;
d(1).cfinv_kin = tmp.cfinv(~d(1).idmec);
d(1).cfinv_mec = tmp.cfinv( d(1).idmec);
end

%% Define initial locking patches
function [d] = InitialLockingPatch(blk,tri,d)
% Coded by Hiroshi Kimura in 2019/2/14
% d(1).idl : Locking, give back-slip  (coupling = 1)
% 
d(1).idl = zeros(blk(1).ntmec,1);
mm1 = 1;
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    nf = size(tri(1).bound(nb1,nb2).clon,2);
    if nf ~= 0
      if blk(1).bound(nb1,nb2).flag2 == 1
        for np = 1:size(blk(1).bound(nb1,nb2).patch,2)
          slipid = inpolygon(tri(1).bound(nb1,nb2).clon,...
                             tri(1).bound(nb1,nb2).clat,...
                             blk(1).bound(nb1,nb2).patch(np).lon,...
                             blk(1).bound(nb1,nb2).patch(np).lat);
          d(1).idl(mm1:mm1+nf-1) = d(1).idl(mm1:mm1+nf-1) | slipid';
        end
        mm1 = mm1 + nf;
      end
    end
  end
end
d(1).idl = logical( d(1).idl);

end

%% Prior function for coupling ratio
function [pdfmc] = prior_mc(mc,minmc,maxmc)
pdfmc = mc<=maxmc & mc>=minmc;
end

%% Prior function for asperity line
function [pdfma] = prior_ma(zu,zd,ZU,ZD,dZ)
% zu, zd: depth (down is +) of asperity lines
% ZU, ZD: limiti depth (down is +) of plate interfaces
% dZ: maximum length (difference) between ZU and ZD (Ortega, 2013)
pdfma = (Heaviside(zu-(ZU-2.*dZ)) - Heaviside(zu- ZD       )).* ...
        (Heaviside(zd- ZU       ) - Heaviside(zd-(ZD+2.*dZ))).* ...
        Heaviside(zd-zu) .* Heaviside(2.*dZ-(ZD-ZU));
end

%% Estimate mechanical coupled area by MCMC (Metropolis-Hasting)
function [cal] = Proceed_MCMC_MH(blk,asp,tri,prm,obs,eul,d,G)
% Test version coded by Hiroshi Kimura in 2019/2/1
% Combined by Hiroshi Kimura in 2019/4/22
% Revised by Hiroshi Kimura in 2019/7/16

% Logging
logfile = fullfile(prm.dirresult,'log.txt');
logfid  = fopen(logfile,'a');
d(1).err = d(1).err./median(d(1).err);
rr = (d(1).obs./d(1).err)'*(d(1).obs./d(1).err);
fprintf('Residual=%9.3f \n',rr);
fprintf(logfid,'Residual=%9.3f \n',rr);

% Center of observation network
alat = mean(obs(1).alat(:));
alon = mean(obs(1).alon(:));

% Heaviside of asperity limit line
Hu = Heaviside(G(1).zulim - G(1).zc);
Hd = Heaviside(G(1).zdlim - G(1).zc);
Hlim = Hd - Hu;

% Initial value
if prm.gpu ~= 99
  precision = 'single'     ;
else
  precision = 'double'     ;
end
rwd         = prm.rwd      ;

% Initial value
mc.int = 1e+1;
mp.int = 1e-10;
mi.int = 1e-10;
ma.int = 1e-5;
la.int = 1e+1;

% Number of unknown parameter
mp.n   = 3.*blk(1).nblock;
mi.n   = 3.*blk(1).nblock;
mc.n   = blk(1).ntkin;
ma.n   = 2.*blk(1).naspline;
ia.n   = blk(1).ntmec;
la.n   = 1;

% Initial std
mp.std = mp.int.*ones(mp.n,1,precision);
mi.std = mi.int.*ones(mi.n,1,precision);
ma.std = ma.int.*ones(ma.n,1,precision);
la.std = la.int.*ones(la.n,1,precision);

% Substitute coupling ratio
mc.old = rand(blk(1).ntkin,1);
% Substitute euler pole vectors
mp.old         = double(blk(1).pole);
mp.old(eul.id) = 0                  ;
mp.old         = mp.old+eul.fixw    ;
% Substitute internal strain tensors
mi.old = 1e-10.*(-0.5+rand(mi.n,1,precision));
mi.old = mi.old.*blk(1).idinter              ;
% Substitute coordinates of up- and down-dip limit
ma.old = zeros(ma.n./2,2);
ma.old(:,1) = blk(1).aline_zu + (blk(1).aline_zd-blk(1).aline_zu).*(0.6 + rand(ma.n./2,1) ./ 5);
ma.old(:,2) = blk(1).aline_zu + (blk(1).aline_zd-blk(1).aline_zu).*(0.1 + rand(ma.n./2,1) ./ 5);
ma.old = reshape(ma.old,ma.n,1);
% Substitute logical value if trimeshes are within asperities or not 
ia.old = zeros(blk(1).ntmec,1);

la.old    = zeros(la.n,1,precision);
res.old   =   inf(   1,1,precision);
% pri.old =   inf(   1,1,precision);

% Scale adjastment of rwd
mcscale  = rwd * 5e-4;
mascale  = rwd * 1e+1;
mpscale  = rwd * 1e-10 .* ones(mp.n,1,precision) .* ~eul.id;
miscale  = rwd * 1e-10;

% Initial chains
cha.mp = zeros(mp.n,prm.kep,precision);
cha.mi = zeros(mi.n,prm.kep,precision);
cha.ma = zeros(ma.n,prm.kep,precision);
cha.mc = zeros(mc.n,prm.kep,precision);
cha.ia = zeros(ia.n,prm.kep,precision);
cha.la = zeros(la.n,prm.kep,precision);

% GPU conversion
if prm.gpu ~= 99
  G(1).tb_mec    = gpuArray(single(full(G(1).tb_mec)));
  G(1).tb_kin    = gpuArray(single(full(G(1).tb_kin)));
  G(1).s         = gpuArray(single(full(G(1).s     )));
  G(1).p         = gpuArray(single(     G(1).p      ));
  G(1).c_kin     = gpuArray(single(     G(1).c_kin  ));
  G(1).c_mec     = gpuArray(single(     G(1).c_mec  ));
  G(1).i         = gpuArray(single(     G(1).i      ));
  d(1).mcid      = gpuArray(single(     d(1).mcid     ));
  d(1).cfinv_mec = gpuArray(single(     d(1).cfinv_mec));
  d(1).cfinv_kin = gpuArray(single(     d(1).cfinv_kin));
end

% MCMC iteration
rt      = 0     ;
count   = 0     ;
burn    = 1     ;
up_mc   = 1     ;
lo_mc   = 0     ;
decrate = 0.9^ 1;
incrate = 0.9^-1;
while not(count == prm.thr)
  rt   = rt+1;
  nacc = 0;tic
  randw = 1;
  % Random value for each parameter
  logu      = log(rand(prm.cha,1,precision));
  rmc       = -randw + (2 * randw) .* rand(mc.n,prm.cha,precision);
  rma       = -randw + (2 * randw) .* rand(ma.n,prm.cha,precision);
  rmp       = -randw + (2 * randw) .* rand(mp.n,prm.cha,precision);
  rmi       = -randw + (2 * randw) .* rand(mi.n,prm.cha,precision);
  rla       = -randw + (2 * randw) .* rand(la.n,prm.cha,precision);

  rmp(         eul.id,:) = 0;
  rmi(~blk(1).idinter,:) = 0;
  for it = 1:prm.cha
    %% Sample section
    %     mctmp            = mc.old+rwd.*mcscale.*rmc(:,it);
    %     id_reject        = mctmp>up_mc | mctmp<lo_mc;
    %     mctmp(id_reject) = mc.old(id_reject);
    %     mc.smp = mctmp                               ;
    %     mp.smp = mp.old + rwd .* mpscale .* rmp(:,it);
    %     mi.smp = mi.old + rwd .* miscale .* rmi(:,it);
    %     la.smp = la.old + rwd .*  la.std .* rla(:,it);
    %     ma.smp = ma.old + rwd .* mascale .* rma(:,it);
    %     %     id_reject = [ma.smp(       1:ma.n/2) < 0 | ma.smp(       1:ma.n/2) > blk(1).aline_zd                                           ;...
    %     %                  ma.smp(ma.n/2+1:   end) < 0 | ma.smp(ma.n/2+1:   end) > blk(1).aline_zd | ma.smp(1:ma.n/2) < ma.smp(ma.n/2+1:end)];
    %     id_reject = [ false(ma.n/2,1)                                                  ;...
    %                  ma.smp(ma.n/2+1:end) < blk(1).aline_zu | ma.smp(ma.n/2+1:end) > blk(1).aline_zd];
    %     ma.smp(id_reject) = ma.old(id_reject);
    %     ma.smp(1:ma.n/2) = max(min(ma.smp(1:ma.n/2),blk(1).aline_zd),ma.smp(ma.n/2+1:end));

    %% Sample section (new method)
    mc.smp = mc.old + rwd .* mcscale .* rmc(:,it);
    mp.smp = mp.old + rwd .* mpscale .* rmp(:,it);
    mi.smp = mi.old + rwd .* miscale .* rmi(:,it);
    la.smp = la.old + rwd .*  la.std .* rla(:,it);
    ma.smp = ma.old + rwd .* mascale .* rma(:,it);
    % Re-sampling coupling ratio
    pdfmc = prior_mc(mc.smp,lo_mc,up_mc);
    while sum(~pdfmc) > 0
      mc.smp(~pdfmc) = mc.old(~pdfmc) + rwd .* mcscale .* (-randw + (2 * randw) .* rand(sum(~pdfmc),1,precision));
      pdfmc = prior_mc(mc.smp,lo_mc,up_mc);
    end
    % Re-sampling asperity line
    pdfma = prior_ma(ma.smp(ma.n/2+1:end),ma.smp(1:ma.n/2),blk(1).aline_zu,blk(1).aline_zd);
    while sum(~pdfma) > 0
      pdfma = repmat(pdfma,2,1);
      ma.smp(~pdfma) = ma.old(~pdfma) + rwd .* mascale .* (-randw + (2 * randw) .* rand(sum(~pdfma),1,precision));
      pdfma = prior_ma(ma.smp(ma.n/2+1:end),ma.smp(1:ma.n/2),blk(1).aline_zu,blk(1).aline_zd);
    end
    
    % Calc gpu memory free capacity
    if prm.gpu ~= 99
      byte1 = whos('G');
      byte2 = whos('mp');
      b = waitGPU(byte1.bytes+byte2.bytes);
    end

    ia.smp = (Heaviside(G(1).zd*ma.smp-G(1).zc) - Heaviside(G(1).zu*ma.smp-G(1).zc)) .* Hlim;
    idl = logical(d(1).maid *  ia.smp);
    idc = logical(d(1).maid * ~ia.smp);
    
    % Calculate back-slip on locked patches.
    bslip              = (G(1).tb_mec * mp.smp) .* d(1).cfinv_mec .* idl;
    
    % Calc inverse Green's function
    %     Gcc        = G(1).s(idc,idc);    % creep -> creep
    %     Gcl        = G(1).s(idc,idl);    % lock  -> creep
    %     bslip(idc) = -Gcc \ (Gcl * bslip(idl));
    bslip(idc) = -G(1).s(idc,idc) \ (G(1).s(idc,idl) * bslip(idl));

    % Due to Rigid motion
    cal.rig = G(1).p * mp.smp;
    % Due to Kinematic coupling
    cal.kin = G(1).c_kin * ((G(1).tb_kin * mp.smp) .* d(1).cfinv_kin .* (d(1).mcid * mc.smp));
    % Due to Mechanical coupling
%     cal.mec = G(1).c_mec * (G(1).E - Gpassive) * bslip;
    cal.mec = G(1).c_mec * bslip;
    % Due to Internal strain
    cal.ine = G(1).i * mi.smp;
    % Zero padding
    if prm.gpu ~= 99
      if isempty(cal.rig); cal.rig = zeros(size(d(1).ind),precision,'gpuArray'); end
      if isempty(cal.kin); cal.kin = zeros(size(d(1).ind),precision,'gpuArray'); end
      if isempty(cal.mec); cal.mec = zeros(size(d(1).ind),precision,'gpuArray'); end
      if isempty(cal.ine); cal.ine = zeros(size(d(1).ind),precision,'gpuArray'); end
    else
      if isempty(cal.rig); cal.rig = zeros(size(d(1).ind)); end
      if isempty(cal.kin); cal.kin = zeros(size(d(1).ind)); end
      if isempty(cal.mec); cal.mec = zeros(size(d(1).ind)); end
      if isempty(cal.ine); cal.ine = zeros(size(d(1).ind)); end
    end        
    % Total velocities
    cal.smp = cal.rig + cal.kin + cal.mec + cal.ine;
    
    if prm.gpu ~= 99
      clear('cal.rig','cal,kin','cal.mec','cal.ine');
    end
    % Calc residual section
    res.smp = sum(((d(1).obs-cal.smp)./d(1).err).^2,1);
    % Mc is better Zero
    %     PRI.SMP = sum(abs(Mc.SMP),1);
    %% MAKE Probably Density Function
    % $$ PDF_{post}=\frac{\frac{1}{\sqrt{2\pi\exp(L)}\times\frac{1}{\sqrt{2\pi}\times\exp{\frac{-Re^{2}}{2}}\exp{\frac{-M^{2}}{2\times\exp{L}}}{\frac{1}{\sqrt{2\pi\exp(L_{old})}\times\frac{1}{\sqrt{2\pi}\times\exp{\frac{-Re^{2}_{old}}{2}}\exp{\frac{-M^{2}_{old}}{2\times\exp{L_{old}}}} $$%%
    %  log(x(x>0));
    %   q1 = logproppdf(x0,y);
    %   q2 = logproppdf(y,x0);
    % This is a generic formula.
    %   rho = (q1+logpdf(y))-(q2+logpdf(x0));
    pdf = -0.5.*...
        ((res.smp+la.smp+exp(-la.smp))...
        -(res.old+la.old+exp(-la.old)));
    %   pdf = -0.5.*(res.smp-res.old);
    acc = pdf > logu(it);
    if acc
      mc.old  =  mc.smp;
      ma.old  =  ma.smp;
      mp.old  =  mp.smp;
      mi.old  =  mi.smp;
      ia.old  =  ia.smp;
      la.old  =  la.smp;
      res.old = res.smp;
      %pri.old = pri.smp;
    end

    % Keep section
    if it > prm.cha - prm.kep
      if prm.gpu ~= 99
        cha.mc(:,it-(prm.cha-prm.kep)) = gather(mc.old);
        cha.ma(:,it-(prm.cha-prm.kep)) = gather(ma.old);
        cha.ia(:,it-(prm.cha-prm.kep)) = gather(ia.old);
        cha.mp(:,it-(prm.cha-prm.kep)) = gather(mp.old);
        cha.mi(:,it-(prm.cha-prm.kep)) = gather(mi.old);
        cha.la(:,it-(prm.cha-prm.kep)) = gather(la.old);
      else
        cha.mc(:,it-(prm.cha-prm.kep)) =        mc.old ;
        cha.ma(:,it-(prm.cha-prm.kep)) =        ma.old ;
        cha.ia(:,it-(prm.cha-prm.kep)) =        ia.old ;
        cha.mp(:,it-(prm.cha-prm.kep)) =        mp.old ;
        cha.mi(:,it-(prm.cha-prm.kep)) =        mi.old ;
        cha.la(:,it-(prm.cha-prm.kep)) =        la.old ;
      end
      if acc; nacc=nacc+1; end
    end
  end
  
% Compress data into int16 format
  CompressData(cha,prm,rt,nacc);
%
  cha.ajr = nacc./prm.cha;
  
% Calc stds from samples
  mc.std = std(cha.mc,1,2);
  ma.std = std(cha.ma,1,2);
  mp.std = std(cha.mp,1,2); 
  mi.std = std(cha.mi,1,2); 
  la.std = std(cha.la,1,2);
  
% Log and display
  fprintf(       't=%3d res=%6.3f accept=%5.1f rwd=%5.2f time=%5.1fsec\n',...
           rt,1-res.old./rr,100*cha.ajr,rwd,toc)
  fprintf(logfid,'t=%3d res=%6.3f accept=%5.1f rwd=%5.2f time=%5.1fsec\n',...
           rt,1-res.old./rr,100*cha.ajr,rwd,toc);
  for bk=1:blk(1).nblock
    [latp,lonp,ang]=xyzp2lla(cha.mp(3.*bk-2,:),cha.mp(3.*bk-1,:),cha.mp(3.*bk,:));
    fprintf(       'pole of block %2i = lat:%7.2f deg. lon:%8.2f deg. ang:%9.2e deg./m.y. \n',...
      bk,mean(latp),mean(lonp),mean(ang));
    fprintf(logfid,'pole of block %2i = lat:%7.2f deg. lon:%8.2f deg. ang:%9.2e deg./m.y. \n',...
      bk,mean(latp),mean(lonp),mean(ang));
  end
  fprintf(       'lamda = %7.2f \n',mean(cha.la));
  fprintf(logfid,'lamda = %7.2f \n',mean(cha.la));

% Adjust random walk distance
  if burn == 0
    if cha.ajr > 0.24
      rwd = rwd * incrate;
    elseif cha.ajr < 0.22
      rwd = rwd * decrate;
    end
    count = count + 1;
  else
    if cha.ajr > 0.24
      rwd = rwd * incrate;
    elseif cha.ajr < 0.22
      rwd = rwd * decrate;
    else
      count = count + 1;
      burn = 0;
    end
  end
  
  cha.smp = cal.smp;
  % Debug-----------
  mpmean = mean(cha.mp,2);
  mcmean = mean(cha.mc,2);
  mamean = mean(cha.ma,2);
  mimean = mean(cha.mi,2);

  iamean = (Heaviside(G(1).zd*mamean-G(1).zc) - Heaviside(G(1).zu*mamean-G(1).zc)) .* Hlim;
  idl = logical(d(1).maid *  iamean);
  idc = logical(d(1).maid * ~iamean);
  
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
% vec.rel = g.c*((g.tb*poltmp).*cf);
  % Debug-----------
  if prm.gpu ~= 99
    ccha.mc  = gather(cha.mc );
    ccha.ma  = gather(cha.ma );
    ccha.mp  = gather(cha.mp );
    ccha.mi  = gather(cha.mi );
    ccha.la  = gather(cha.la );
    ccha.smp = gather(cha.smp);
    MakeFigures(ccha,blk,obs,rt,gather(lo_mc),gather(up_mc),vec,bslip,bslipl,mimean)
  else
    MakeFigures( cha,blk,obs,rt,       lo_mc ,       up_mc ,vec,bslip,bslipl,mimean)
  end
  if rt > prm.itr; break; end
end

% Display
cha.res = res.smp;
fprintf(       'RMS=: %8.3f\n',cha.res)
fprintf(logfid,'RMS=: %8.3f\n',cha.res);
fprintf(       '=== Finished MCMC part ===\n')
fprintf(logfid,'=== finished MCMC part ===\n');
fclose(logfid);
end

%% Estimate mechanical coupled area by MCMC (Replica exchange method)
function [cal] = Proceed_MCMC_RE(blk,asp,tri,prm,obs,eul,d,G)
% Test version coded by Hiroshi Kimura in 2019/2/1
% Combined by Hiroshi Kimura in 2019/4/22
% Revised by Hiroshi Kimura in 2019/7/16

% Logging
logfile = fullfile(prm.dirresult,'log.txt');
logfid  = fopen(logfile,'a');
d(1).invms = 1 ./ ( d(1).obs'*d(1).obs ./ (size(d(1).obs,1)/3) );
d(1).sigma = sqrt(d(1).err.^2 + ones(size(d(1).err)).^2);
rr = (d(1).obs./d(1).err)'*(d(1).obs./d(1).err);
rr_inv =  1 / rr;
fprintf('Residual=%9.3f \n',rr);
fprintf(logfid,'Residual=%9.3f \n',rr);

% Center of observation network
alat = mean(obs(1).alat(:));
alon = mean(obs(1).alon(:));

% Heaviside of asperity limit line
Hu = repmat(Heaviside(G(1).zulim - G(1).zc),1,prm.nrep);
Hd = repmat(Heaviside(G(1).zdlim - G(1).zc),1,prm.nrep);
Hlim = Hd - Hu;

% Initial value
if prm.gpu ~= 99
  precision = 'single'     ;
else
  precision = 'double'     ;
end
rwd         = prm.rwd      ;

% Inverse temperatures
Tmax  = prm.tmax;
dT    = Tmax^(1/(prm.nrep-1));
T_inv = 1 ./ (dT.^[0:prm.nrep-1]);

% Initial value
mc.int = 1e+1;
mp.int = 1e-10;
mi.int = 1e-10;
ma.int = 1e-5;
la.int = 1e+1;

% Number of unknown parameter
mp.n   = 3.*blk(1).nblock;
mi.n   = 3.*blk(1).nblock;
mc.n   = blk(1).ntkin;
ma.n   = 2.*blk(1).naspline;
ia.n   = blk(1).ntmec;
la.n   = 1;

% Initial std
mp.std = mp.int.*ones(mp.n,prm.nrep,precision);
mi.std = mi.int.*ones(mi.n,prm.nrep,precision);
ma.std = ma.int.*ones(ma.n,prm.nrep,precision);
la.std = la.int.*ones(la.n,prm.nrep,precision);

% Substitute coupling ratio
mc.old = rand(blk(1).ntkin,prm.nrep);
% Substitute euler pole vectors
mp.old           = repmat(double(blk(1).pole),1,prm.nrep);
mp.old(eul.id,:) = 0                                     ;
mp.old           = mp.old+repmat(eul.fixw,1,prm.nrep)    ;
% Substitute internal strain tensors
mi.old = 1e-10.*(-0.5+rand(mi.n,prm.nrep,precision));
mi.old = mi.old.*repmat(blk(1).idinter,1,prm.nrep)  ;
% Substitute coordinates of up- and down-dip limit
ma.old = zeros(ma.n,prm.nrep);
ma.old(       1:ma.n/2,:) = blk(1).aline_zu + repmat(blk(1).aline_dz,1,prm.nrep).*(0.9 + rand(ma.n./2,prm.nrep) ./ 5);
ma.old(ma.n/2+1:   end,:) = blk(1).aline_zu + repmat(blk(1).aline_dz,1,prm.nrep).*(0.1 + rand(ma.n./2,prm.nrep) ./ 5);
% Substitute logical value if trimeshes are within asperities or not 
ia.old = zeros(blk(1).ntmec,prm.nrep);

la.old    =  ones(la.n,prm.nrep,precision);
res.old   =   inf(   1,prm.nrep,precision);
res.oldl  =   inf(   1,prm.nrep,precision);
% pri.old   =   inf(   1,prm.nrep,precision);

% Scale adjastment of rwd
mcscale  = rwd * 5e-4;
mascale  = rwd * 3e-1 .* (1e-2 .* repmat(blk(1).aline_dz,2,prm.nrep));
mpscale  = rwd * 1e-10 .* repmat(ones(mp.n,1,precision).*~eul.id,1,prm.nrep);
miscale  = rwd * 1e-10;
lascale  = rwd * 1e+0;

% Initial chains
cha.mp = zeros(mp.n,prm.kep,prm.nrep,precision);
cha.mi = zeros(mi.n,prm.kep,prm.nrep,precision);
cha.ma = zeros(ma.n,prm.kep,prm.nrep,precision);
cha.mc = zeros(mc.n,prm.kep,prm.nrep,precision);
cha.ia = zeros(ia.n,prm.kep,prm.nrep,precision);
cha.la = zeros(la.n,prm.kep,prm.nrep,precision);

% GPU conversion
if prm.gpu ~= 99
  G(1).tb_mec    = gpuArray(single(full(G(1).tb_mec)));
  G(1).tb_kin    = gpuArray(single(full(G(1).tb_kin)));
  G(1).s         = gpuArray(single(full(G(1).s     )));
  G(1).p         = gpuArray(single(     G(1).p      ));
  G(1).c_kin     = gpuArray(single(     G(1).c_kin  ));
  G(1).c_mec     = gpuArray(single(     G(1).c_mec  ));
  G(1).i         = gpuArray(single(     G(1).i      ));
  d(1).mcid      = gpuArray(single(     d(1).mcid     ));
  d(1).cfinv_mec = gpuArray(single(     d(1).cfinv_mec));
  d(1).cfinv_kin = gpuArray(single(     d(1).cfinv_kin));
  rr_inv         = gpuArray(single( rr_inv));
  res.old        = gpuArray(single(res.old));
  res.oldl       = gpuArray(single(res.oldl));
end

% MCMC iteration
rt      = 0     ;
count   = 0     ;
burn    = 1     ;
up_mc   = 1     ;
lo_mc   = 0     ;
decrate = 0.9^ 1;
incrate = 0.9^-1;
while not(count == prm.thr)
  rt   = rt+1;
  nacc = zeros(1,prm.nrep);tic
  randw = 1;
  excount = zeros(1,prm.nrep);
  % Random value for each parameter
  logu      = log(rand(prm.cha,prm.nrep  ,precision));
  loge      = log(rand(prm.cha,prm.nrep-1,precision));
  rmc       = -randw + (2 * randw) .* rand(mc.n,prm.nrep,prm.cha,precision);
  rma       = -randw + (2 * randw) .* rand(ma.n,prm.nrep,prm.cha,precision);
  rmp       = -randw + (2 * randw) .* rand(mp.n,prm.nrep,prm.cha,precision);
  rmi       = -randw + (2 * randw) .* rand(mi.n,prm.nrep,prm.cha,precision);
  rla       = -randw + (2 * randw) .* rand(la.n,prm.nrep,prm.cha,precision);
  rex       = randi(prm.nrep-1,1,prm.cha);
  
  rmp(         eul.id,:,:) = 0;
  rmi(~blk(1).idinter,:,:) = 0;
  for it = 1:prm.cha
    % Sample section
    mc.smp = mc.old + rwd .* mcscale .* rmc(:,:,it);
    mp.smp = mp.old + rwd .* mpscale .* rmp(:,:,it);
    mi.smp = mi.old + rwd .* miscale .* rmi(:,:,it);
    la.smp = la.old + rwd .* lascale .* rla(:,:,it);
    ma.smp = ma.old + rwd .* mascale .* rma(:,:,it);
    % Re-sampling coupling ratio
    pdfmc = prior_mc(mc.smp,lo_mc,up_mc);
    while sum(sum(~pdfmc)) > 0
      mc.smp(~pdfmc) = mc.old(~pdfmc) + rwd .* mcscale .* (-randw + (2 * randw) .* rand(sum(sum(~pdfmc)),1,precision));
      pdfmc = prior_mc(mc.smp,lo_mc,up_mc);
    end
    % Re-sampling asperity line
    pdfma = prior_ma(ma.smp(ma.n/2+1:end,:),ma.smp(1:ma.n/2,:),blk(1).aline_zu,blk(1).aline_zd,blk(1).aline_dz);
    while sum(sum(~pdfma)) > 0
      pdfma = repmat(pdfma,2,1);
      ma.smp(~pdfma) = ma.old(~pdfma) + rwd .* mascale(~pdfma) .* (-randw + (2 * randw) .* rand(sum(sum(~pdfma)),1,precision));
      pdfma = prior_ma(ma.smp(ma.n/2+1:end,:),ma.smp(1:ma.n/2,:),blk(1).aline_zu,blk(1).aline_zd,blk(1).aline_dz);
    end
    % Re-sampling lamda
    %     pdfla = prior_mc(la.smp,0,Inf);
    %     while sum(sum(~pdfla)) > 0
    %       la.smp(~pdfla) = la.old(~pdfla) + rwd .* lascale .* (-randw + (2 * randw) .* rand(1,sum(sum(~pdfla)),precision));
    %       pdfla = prior_mc(la.smp,0,Inf);
    %     end
    % Calc gpu memory free capacity
    if prm.gpu ~= 99
      byte1 = whos('G');
      byte2 = whos('mp');
      b = waitGPU(byte1.bytes+byte2.bytes);
    end

    ia.smp = (Heaviside(G(1).zd*ma.smp-G(1).zc) - Heaviside(G(1).zu*ma.smp-G(1).zc)) .* Hlim;
    idl    = logical(d(1).maid *  ia.smp);
    idc    = logical(d(1).maid * ~ia.smp);
    
    % Calculate back-slip on locked patches.
    bslip              = (G(1).tb_mec * mp.smp) .* d(1).cfinv_mec .* idl;
    
    % Calc inverse Green's function
    %     Gcc        = G(1).s(idc,idc);    % creep -> creep
    %     Gcl        = G(1).s(idc,idl);    % lock  -> creep
    %     bslip(idc) = -Gcc \ (Gcl * bslip(idl));
    for nrep=1:prm.nrep
      bslip(idc(:,nrep),nrep) = -G(1).s(idc(:,nrep),idc(:,nrep)) \ (G(1).s(idc(:,nrep),idl(:,nrep)) * bslip(idl(:,nrep),nrep));
    end
    
    % Due to Rigid motion
    cal.rig = G(1).p * mp.smp;
    % Due to Kinematic coupling
    cal.kin = G(1).c_kin * ((G(1).tb_kin * mp.smp) .* d(1).cfinv_kin .* (d(1).mcid * mc.smp));
    % Due to Mechanical coupling
%     cal.mec = G(1).c_mec * (G(1).E - Gpassive) * bslip;
    cal.mec = G(1).c_mec * bslip;
    % Due to Internal strain
    cal.ine = G(1).i * mi.smp;
    % Zero padding
    if prm.gpu ~= 99
      if isempty(cal.rig); cal.rig = zeros(size(d(1).ind,1),prm.nrep,precision,'gpuArray'); end
      if isempty(cal.kin); cal.kin = zeros(size(d(1).ind,1),prm.nrep,precision,'gpuArray'); end
      if isempty(cal.mec); cal.mec = zeros(size(d(1).ind,1),prm.nrep,precision,'gpuArray'); end
      if isempty(cal.ine); cal.ine = zeros(size(d(1).ind,1),prm.nrep,precision,'gpuArray'); end
    else
      if isempty(cal.rig); cal.rig = zeros(size(d(1).ind,1),prm.nrep); end
      if isempty(cal.kin); cal.kin = zeros(size(d(1).ind,1),prm.nrep); end
      if isempty(cal.mec); cal.mec = zeros(size(d(1).ind,1),prm.nrep); end
      if isempty(cal.ine); cal.ine = zeros(size(d(1).ind,1),prm.nrep); end
    end        
    % Total velocities
    cal.smp = cal.rig + cal.kin + cal.mec + cal.ine;
    
    if prm.gpu ~= 99
      clear('cal.rig','cal,kin','cal.mec','cal.ine');
    end
    % Calc residual section
    res.smp = sum(((d(1).obs-cal.smp)./d(1).err).^2,1);
    %     res.smpl= sum(((d(1).obs-cal.smp)./d(1).err).^2,1) ./ (la.smp).^2;  % considering model error
    res.smpl= sum(abs(d(1).obs-cal.smp)./(d(1).sigma),1);  % L-1 norm (Ortega, 2013)
    % Mc is better Zero
    %% MAKE Probably Density Function
    %     pdf = -0.5 .* (res.smp-res.old) .* rr_inv .* T_inv;   % normalized by obs errors
    %     pdf = (-2 .* (log(la.smp) - log(la.old)) -0.5 .* (res.smpl - res.oldl)) .* T_inv;  % consider model uncertainties
    pdf = -1 .* (res.smpl - res.oldl) .* T_inv;   % PDF of L1-norm (Ortega, 2013)
    
    % Accept 
    acc = pdf > logu(it,:);
    mc.old( :,acc) = mc.smp( :,acc);
    ma.old( :,acc) = ma.smp( :,acc);
    mp.old( :,acc) = mp.smp( :,acc);
    mi.old( :,acc) = mi.smp( :,acc);
    ia.old( :,acc) = ia.smp( :,acc);
    la.old( :,acc) = la.smp( :,acc);
    res.old(:,acc) = res.smp(:,acc);
    res.oldl(:,acc)= res.smpl(:,acc);
    %     pri.old(:,acc)  = pri.smp(:,acc) ;

    % Exchange Replicas
    if mod(it,prm.efrq) == 0
      %       r = -0.5 .* (res.old(rex(it)+1)-res.old(rex(it))) * rr_inv * (T_inv(rex(it))-T_inv(rex(it)+1));
      %       r = -2 .* (T_inv(rex(it))-T_inv(rex(it)+1)) .* (log(la.smp(rex(it)+1)) - log(la.smp(rex(it)))) -0.5 .* (res.oldl(rex(it)+1) - res.oldl(rex(it))) .* (T_inv(rex(it))-T_inv(rex(it)+1));
      %       r = -(res.oldl(rex(it)+1)-res.oldl(rex(it))) .* (T_inv(rex(it))-T_inv(rex(it)+1));
      %       if r > loge(it)
      %         mc.old(:,[rex(it),rex(it)+1]) = fliplr(mc.old(:,[rex(it),rex(it)+1]));
      %         ma.old(:,[rex(it),rex(it)+1]) = fliplr(ma.old(:,[rex(it),rex(it)+1]));
      %         mp.old(:,[rex(it),rex(it)+1]) = fliplr(mp.old(:,[rex(it),rex(it)+1]));
      %         mi.old(:,[rex(it),rex(it)+1]) = fliplr(mi.old(:,[rex(it),rex(it)+1]));
      %         ia.old(:,[rex(it),rex(it)+1]) = fliplr(ia.old(:,[rex(it),rex(it)+1]));
      %         la.old(:,[rex(it),rex(it)+1]) = fliplr(la.old(:,[rex(it),rex(it)+1]));
      %         res.old(:,[rex(it),rex(it)+1]) = fliplr(res.old(:,[rex(it),rex(it)+1]));
      %         res.oldl(:,[rex(it),rex(it)+1]) = fliplr(res.oldl(:,[rex(it),rex(it)+1]));
      %         excount(nrep) = excount(nrep) + 1;
      %       end
      for nrep = 1:prm.nrep-1
        %         r = -0.5 .* (res.old(nrep+1)-res.old(nrep)) * rr_inv * (T_inv(nrep)-T_inv(nrep+1));
        %         r = -2 .* (T_inv(rex(it))-T_inv(rex(it)+1)) .* (log(la.smp(rex(it)+1)) - log(la.smp(rex(it)))) -0.5 .* (res.oldl(rex(it)+1) - res.oldl(rex(it))) .* (T_inv(rex(it))-T_inv(rex(it)+1));
        r = -(res.oldl(nrep+1)-res.oldl(nrep)) .* (T_inv(nrep)-T_inv(nrep+1));
        if r > loge(it,nrep)
          mc.old(:,[nrep,nrep+1]) = fliplr(mc.old(:,[nrep,nrep+1]));
          ma.old(:,[nrep,nrep+1]) = fliplr(ma.old(:,[nrep,nrep+1]));
          mp.old(:,[nrep,nrep+1]) = fliplr(mp.old(:,[nrep,nrep+1]));
          mi.old(:,[nrep,nrep+1]) = fliplr(mi.old(:,[nrep,nrep+1]));
          ia.old(:,[nrep,nrep+1]) = fliplr(ia.old(:,[nrep,nrep+1]));
          la.old(:,[nrep,nrep+1]) = fliplr(la.old(:,[nrep,nrep+1]));
          res.old(:,[nrep,nrep+1]) = fliplr(res.old(:,[nrep,nrep+1]));
          res.oldl(:,[nrep,nrep+1])= fliplr(res.oldl(:,[nrep,nrep+1]));
          excount(nrep) = excount(nrep) + 1;
        end
      end
    end

    % Keep section
    if it > prm.cha - prm.kep
      if prm.gpu ~= 99
        cha.mc(:,it-(prm.cha-prm.kep),:) = gather(mc.old);
        cha.ma(:,it-(prm.cha-prm.kep),:) = gather(ma.old);
        cha.ia(:,it-(prm.cha-prm.kep),:) = gather(ia.old);
        cha.mp(:,it-(prm.cha-prm.kep),:) = gather(mp.old);
        cha.mi(:,it-(prm.cha-prm.kep),:) = gather(mi.old);
        cha.la(:,it-(prm.cha-prm.kep),:) = gather(la.old);
      else
        cha.mc(:,it-(prm.cha-prm.kep),:) =        mc.old ;
        cha.ma(:,it-(prm.cha-prm.kep),:) =        ma.old ;
        cha.ia(:,it-(prm.cha-prm.kep),:) =        ia.old ;
        cha.mp(:,it-(prm.cha-prm.kep),:) =        mp.old ;
        cha.mi(:,it-(prm.cha-prm.kep),:) =        mi.old ;
        cha.la(:,it-(prm.cha-prm.kep),:) =        la.old ;
      end
      nacc(acc) = nacc(acc) + 1;
    end
    
  end
  
  % Compress data into int16 format
  CompressData(cha,prm,rt,nacc);
  %
  cha.ajr = nacc./prm.cha;
  
  % Calc stds from samples
  mc.std = std(cha.mc,1,2);
  ma.std = std(cha.ma,1,2);
  mp.std = std(cha.mp,1,2); 
  mi.std = std(cha.mi,1,2); 
  la.std = std(cha.la,1,2);
  
  % Log and display
  fprintf(       '\nt=%3d rwd=%5.2f time=%5.1sec\n',rt,rwd,toc)
  fprintf(logfid,'\nt=%3d rwd=%5.2f time=%5.1sec\n',rt,rwd,toc);
  for nrep = 1:prm.nrep
      fprintf(       'N=%2i res=%6.3f accept=%5.1f%% lamda=%7.2f exratio=%5.2f%%\n',...
           nrep,1-res.old(nrep)./rr,100*cha.ajr(nrep),mean(cha.la(:,:,nrep)),100*excount(nrep)/prm.cha)
      fprintf(logfid,'N=%2i res=%6.3f accept=%5.1f%% lamda=%7.2f exratio=%5.2f%%\n',...
           nrep,1-res.old(nrep)./rr,100*cha.ajr(nrep),mean(cha.la(:,:,nrep)),100*excount(nrep)/prm.cha);
  end
  for bk=1:blk(1).nblock
    [latp,lonp,ang]=xyzp2lla(cha.mp(3.*bk-2,:),cha.mp(3.*bk-1,:),cha.mp(3.*bk,:));
    fprintf(       'pole of block %2i = lat:%7.2f deg. lon:%8.2f deg. ang:%9.2e deg./m.y. \n',...
      bk,mean(latp),mean(lonp),mean(ang));
    fprintf(logfid,'pole of block %2i = lat:%7.2f deg. lon:%8.2f deg. ang:%9.2e deg./m.y. \n',...
      bk,mean(latp),mean(lonp),mean(ang));
  end

  % Adjust random walk distance
  if burn == 0
    if cha.ajr(1) > 0.24
      rwd = rwd * incrate;
    elseif cha.ajr(1) < 0.22
      rwd = rwd * decrate;
    end
    count = count + 1;
  else
    if cha.ajr(1) > 0.24
      rwd = rwd * incrate;
    elseif cha.ajr(1) < 0.22
      rwd = rwd * decrate;
    else
      count = count + 1;
      burn = 0;
    end
  end
  
  cha.smp = gather(cal.smp);
  % Debug-----------
  mpmean = permute(mean(cha.mp,2),[1 3 2]);
  mcmean = permute(mean(cha.mc,2),[1 3 2]);
  mamean = permute(mean(cha.ma,2),[1 3 2]);
  mimean = permute(mean(cha.mi,2),[1 3 2]);

  iamean = (Heaviside(G(1).zd*mamean-G(1).zc) - Heaviside(G(1).zu*mamean-G(1).zc)) .* Hlim;
  idl = logical(d(1).maid *  iamean);
  idc = logical(d(1).maid * ~iamean);
  
  % Calculate back-slip on locked and creeping patches.
  bslip      = (G(1).tb_mec * mpmean) .* d(1).cfinv_mec .* idl;
  bslipl     = bslip;
  for nrep=1:prm.nrep
    bslip(idc(:,nrep),nrep) = -G(1).s(idc(:,nrep),idc(:,nrep)) \ (G(1).s(idc(:,nrep),idl(:,nrep)) * bslip(idl(:,nrep),nrep));
  end

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
  % Debug-----------
  if prm.gpu ~= 99
    ccha.mc  = gather(cha.mc );
    ccha.ma  = gather(cha.ma );
    ccha.ia  = gather(cha.ia );
    ccha.mp  = gather(cha.mp );
    ccha.mi  = gather(cha.mi );
    ccha.la  = gather(cha.la );
    ccha.smp = gather(cha.smp);
    MakeFigures(ccha,blk,obs,rt,gather(lo_mc),gather(up_mc),vec,bslip,bslipl,mimean)
  else
    MakeFigures( cha,blk,obs,rt,       lo_mc ,       up_mc ,vec,bslip,bslipl,mimean)
  end
  if rt > prm.itr; break; end
end

% Display
cha.res = res.smp;
fprintf(       'RMS=: %8.3f\n',cha.res(1))
fprintf(logfid,'RMS=: %8.3f\n',cha.res(1));
fprintf(       '=== Finished MCMC part ===\n')
fprintf(logfid,'=== finished MCMC part ===\n');
fclose(logfid);
end

%% Estimate passive slip distribution by forward simulation
function [cal] = CalcPassiveSlip(blk,asp,tri,prm,obs,eul,d,G)
precision = 'single';
% initial parameters
mc.int = 1e-2 ;
mp.int = 1e-10;
mi.int = 1e-10;
mc.n   = blk(1).ntkin;
mp.n   = 3.*blk(1).nblock;
mi.n   = 3.*blk(1).nblock;
mc.old = zeros(mc.n,1);
% substitute euler pole vectors
mp.old         = double(blk(1).pole);
mp.old(eul.id) = 0                  ;
mp.old         = mp.old+eul.fixw    ;
% substitute internal strain tensors
mi.old = 1e-10.*(-0.5+rand(mi.n,1,precision));
mi.old = mi.old.*blk(1).idinter              ;

% substitude index of locking patches
idl = logical(d(1).maid *  d(1).idl);
idc = logical(d(1).maid * ~d(1).idl);

% Calculate back-slip on locked patches.
bslip              = (G(1).tb_mec * mp.old) .* d(1).cfinv_mec .* idl;
cal.aslip = bslip;

% Calc inverse Green's function
%     Gcc        = G(1).s(idc,idc);    % creep -> creep
%     Gcl        = G(1).s(idc,idl);    % lock  -> creep
%     bslip(idc) = -Gcc \ (Gcl * bslip(idl));
bslip(idc) = -G(1).s(idc,idc) \ (G(1).s(idc,idl) * bslip(idl));
cal.bslip = bslip;

% Due to Rigid motion
cal.rig = G(1).p * mp.old;
% Due to Kinematic coupling
cal.kin = G(1).c_kin * ((G(1).tb_kin * mp.old) .* d(1).cfinv_kin .* (d(1).mcid * mc.old));
% Due to Mechanical coupling
%     cal.mec = G(1).c_mec * (G(1).E - Gpassive) * bslip;
cal.mec = G(1).c_mec * bslip;
% Due to Internal strain
cal.ine = G(1).i * mi.old;
% Zero padding
if isempty(cal.rig); cal.rig = zeros(size(d(1).ind)); end
if isempty(cal.kin); cal.kin = zeros(size(d(1).ind)); end
if isempty(cal.mec); cal.mec = zeros(size(d(1).ind)); end
if isempty(cal.ine); cal.ine = zeros(size(d(1).ind)); end
% Total velocities
cal.smp = cal.rig + cal.kin + cal.mec + cal.ine;

MakeFigs(blk,cal,bslip,obs)

end

%% Compress CHA sampling
function CompressData(ccha,prm,itr,nacc)
% Compressing CHA sampled parameter to int
% sfactor = 2^8 ;  % int8
sfactor = 2^16;  % int16
% 
ccha.mc =  single(ccha.mc);
ccha.ma =  single(ccha.ma);
ccha.ia = logical(ccha.ia);
ccha.mp =  single(ccha.mp);
ccha.mi =  single(ccha.mi);
% if prm.gpu==99&&gpudevicecount==0
if prm.gpu == 99
  meanmc = mean(ccha.mc,2);
  meanma = mean(ccha.ma,2);
  meania = mean(ccha.ia,2);
  meanmp = mean(ccha.mp,2);
  meanmi = mean(ccha.mi,2);
  for nrep = 1:prm.nrep
    covmc(:,:,nrep) = cov(ccha.mc(:,:,nrep)');
    covma(:,:,nrep) = cov(ccha.ma(:,:,nrep)');
    covia(:,:,nrep) = cov(ccha.ia(:,:,nrep)');
    covmp(:,:,nrep) = cov(ccha.mp(:,:,nrep)');
    covmi(:,:,nrep) = cov(ccha.mi(:,:,nrep)');  
  end
else
  gcha.mc = gpuArray(ccha.mc);
  gcha.ma = gpuArray(ccha.ma);
  gcha.ia = gpuArray(ccha.ia);
  gcha.mp = gpuArray(ccha.mp);
  gcha.mi = gpuArray(ccha.mi);
  meanmc  =  mean(gcha.mc,2);
  meanma  =  mean(gcha.ma,2);
  meania  =  mean(gcha.ia,2);
  meanmp  =  mean(gcha.mp,2);
  meanmi  =  mean(gcha.mi,2);
  for nrep = 1:prm.nrep
    covmc(:,:,nrep) = cov(gcha.mc(:,:,nrep)');
    covma(:,:,nrep) = cov(gcha.ma(:,:,nrep)');
    covia(:,:,nrep) = cov(gcha.ia(:,:,nrep)');
    covmp(:,:,nrep) = cov(gcha.mp(:,:,nrep)');
    covmi(:,:,nrep) = cov(gcha.mi(:,:,nrep)');
  end
  meanmc = gather(meanmc);
  meanma = gather(meanma);
  meania = gather(meania);
  meanmp = gather(meanmp);
  meanmi = gather(meanmi);
  covmc  =  gather(covmc);
  covma  =  gather(covma);
  covia  =  gather(covia);
  covmp  =  gather(covmp);
  covmi  =  gather(covmi);
end
% 
mcmax = max(ccha.mc,[],2);
mcmin = min(ccha.mc,[],2);
mamax = max(ccha.ma,[],2);
mamin = min(ccha.ma,[],2);
iamax = max(ccha.ia,[],2);
iamin = min(ccha.ia,[],2);
mpmax = max(ccha.mp,[],2);
mpmin = min(ccha.mp,[],2);
mimax = max(ccha.mi,[],2);
mimin = min(ccha.mi,[],2);
% 
mcscale = 1./(mcmax-mcmin);
mcbase  = bsxfun(@minus,bsxfun(@times,bsxfun(@minus,ccha.mc,mcmin),mcscale.*(sfactor-1)),sfactor/2);
% mcint = int8(mcbase);
mcint   = int16(mcbase);

mascale = 1./(mamax-mamin);
mabase  = bsxfun(@minus,bsxfun(@times,bsxfun(@minus,ccha.ma,mamin),mascale.*(sfactor-1)),sfactor/2);
% maint = int8(mabase);
maint   = int16(mabase);

mpscale = 1./(mpmax-mpmin);
mpbase  = bsxfun(@minus,bsxfun(@times,bsxfun(@minus,ccha.mp,mpmin),mpscale.*(sfactor-1)),sfactor/2);
% mpint = int8(mpbase);
mpint   = int16(mpbase);

miscale = 1./(mimax-mimin);
mibase  = bsxfun(@minus,bsxfun(@times,bsxfun(@minus,ccha.mi,mimin),miscale.*(sfactor-1)),sfactor/2);
% miint = int8(mibase);
miint   = int16(mibase);

% mc
for ii = 1:size(mcint,1)
  cha.mccompress.nflt(ii).mcscale = mcscale(ii,:,:);
  cha.mccompress.nflt(ii).mcmax   =   mcmax(ii,:,:);
  cha.mccompress.nflt(ii).mcmin   =   mcmin(ii,:,:);
end
cha.mccompress.covmc   =         covmc;
cha.mccompress.meanmc  =        meanmc;
% cha.mccompress.smpmc =  int8(mcbase);
cha.mccompress.smpmc   = int16(mcbase);

% ma
for ii = 1:size(maint,1)
  cha.macompress.nasp(ii).mascale = mascale(ii,:,:);
  cha.macompress.nasp(ii).mamax   =   mamax(ii,:,:);
  cha.macompress.nasp(ii).mamin   =   mamin(ii,:,:);
end
cha.macompress.covma   =         covma;
cha.macompress.meanma  =        meanma;
% cha.macompress.smpma =  int8(mabase);
cha.macompress.smpma   = int16(mabase);

% ia
cha.iacompress.covia  =            covia;
cha.iacompress.meania =           meania;
cha.iacompress.smpia  = logical(ccha.ia);

% mp
for ii = 1:size(mpint,1)
  cha.mpcompress.npol(ii).mpscale = mpscale(ii,:,:);
  cha.mpcompress.npol(ii).mpmax   =   mpmax(ii,:,:);
  cha.mpcompress.npol(ii).mpmin   =   mpmin(ii,:,:);
end
cha.mpcompress.covmp   =         covmp;
cha.mpcompress.meanmp  =        meanmp;
% cha.mpcompress.smpmp =  int8(mpbase);
cha.mpcompress.smpmp   = int16(mpbase);

% mi
for ii = 1:size(miint,1)
  cha.micompress.nine(ii).miscale = miscale(ii,:,:);
  cha.micompress.nine(ii).mimax   =   mimax(ii,:,:);
  cha.micompress.nine(ii).mimin   =   mimin(ii,:,:);
end
cha.micompress.covmi   =         covmi;
cha.micompress.meanmi  =        meanmi;
% cha.micompress.smpmi =  int8(mibase);
cha.micompress.smpmi   = int16(mibase);

cha.ajr = nacc./prm.cha;
save(fullfile(prm.dirresult,['cha_test',num2str(itr,'%03i')]),'cha','-v7.3');
%
end

%% Define random walk lines
function [blk,asp] = DefRandomWalkLine(blk,prm,obs)

asp  = [];
alat = mean(obs(1).alat(:));
alon = mean(obs(1).alon(:));
blk(1).naspline   =  0;
blk(1).aline_zu   = [];
blk(1).aline_zd   = [];
blk(1).aline_dz   = [];
blk(1).aline_lonu = [];
blk(1).aline_lond = [];
blk(1).aline_latu = [];
blk(1).aline_latd = [];

Nint = 1;  % default

[tricx,tricy] = PLTXY(blk(1).clat,blk(1).clon,alat,alon);
% 
mm3 = 1;
mm1 = 1;
np = 1;
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    nf = size(blk(1).bound(nb1,nb2).blon,1);
    if nf ~= 0
      if blk(1).bound(nb1,nb2).flag2 == 1
        if size(blk(1).bound(nb1,nb2).udshare,1) == 0
          rwlfile = fullfile(prm.dirblock_interface,['udlineb_',num2str(nb1),'_',num2str(nb2),'.txt']);
        else
          rwlfile = fullfile(prm.dirblock_interface,['udlineb_',num2str(blk(1).bound(nb1,nb2).udshare(1)),'_',num2str(blk(1).bound(nb1,nb2).udshare(2)),'.txt']);
        end
        fid     = fopen(rwlfile,'r');
        if fid >= 0
          asp(np).nb1 = nb1;
          asp(np).nb2 = nb2;
          if size(blk(1).bound(nb1,nb2).udinterp,1) ~= 0
            nint = blk(1).bound(nb1,nb2).udinterp;
          else
            nint = Nint;
          end
          tmp = fscanf(fid,'%f %f %f %f %f %f \n',[6, Inf]);
          nasp = size(tmp,2);
          blk(1).bound(nb1,nb2).asp_lond = tmp(1,:)';
          blk(1).bound(nb1,nb2).asp_latd = tmp(2,:)';
          blk(1).bound(nb1,nb2).asp_depd = tmp(3,:)';
          blk(1).bound(nb1,nb2).asp_lonu = tmp(4,:)';
          blk(1).bound(nb1,nb2).asp_latu = tmp(5,:)';
          blk(1).bound(nb1,nb2).asp_depu = tmp(6,:)';
          blk(1).bound(nb1,nb2).asp_ddep = abs(blk(1).bound(nb1,nb2).asp_depd - blk(1).bound(nb1,nb2).asp_depu);
          [xd,yd] = PLTXY(blk(1).bound(nb1,nb2).asp_latd,blk(1).bound(nb1,nb2).asp_lond,alat,alon);
          [xu,yu] = PLTXY(blk(1).bound(nb1,nb2).asp_latu,blk(1).bound(nb1,nb2).asp_lonu,alat,alon);
          blk(1).bound(nb1,nb2).asp_xd =                             xd;
          blk(1).bound(nb1,nb2).asp_yd =                             yd;
          blk(1).bound(nb1,nb2).asp_zd = blk(1).bound(nb1,nb2).asp_depd;
          blk(1).bound(nb1,nb2).asp_xu =                             xu;
          blk(1).bound(nb1,nb2).asp_yu =                             yu;
          blk(1).bound(nb1,nb2).asp_zu = blk(1).bound(nb1,nb2).asp_depu;
          % Interpolation of up- and down-dip point
          blk(1).bound(nb1,nb2).asp_xd_interp = linspace2(blk(1).bound(nb1,nb2).asp_xd, nint);
          blk(1).bound(nb1,nb2).asp_yd_interp = linspace2(blk(1).bound(nb1,nb2).asp_yd, nint);
          blk(1).bound(nb1,nb2).asp_zd_interp = linspace2(blk(1).bound(nb1,nb2).asp_zd, nint);
          blk(1).bound(nb1,nb2).asp_xu_interp = linspace2(blk(1).bound(nb1,nb2).asp_xu, nint);
          blk(1).bound(nb1,nb2).asp_yu_interp = linspace2(blk(1).bound(nb1,nb2).asp_yu, nint);
          blk(1).bound(nb1,nb2).asp_zu_interp = linspace2(blk(1).bound(nb1,nb2).asp_zu, nint);
          np = np + 1;
          if size(blk(1).bound(nb1,nb2).udshare,1) == 0
            blk(1).naspline  =  blk(1).naspline + nasp;
            blk(1).aline_zu   = [blk(1).aline_zu  ; blk(1).bound(nb1,nb2).asp_depu];
            blk(1).aline_zd   = [blk(1).aline_zd  ; blk(1).bound(nb1,nb2).asp_depd];
            blk(1).aline_dz   = [blk(1).aline_dz  ; blk(1).bound(nb1,nb2).asp_ddep];
            blk(1).aline_lonu = [blk(1).aline_lonu; blk(1).bound(nb1,nb2).asp_lonu];
            blk(1).aline_lond = [blk(1).aline_lond; blk(1).bound(nb1,nb2).asp_lond];
            blk(1).aline_latu = [blk(1).aline_latu; blk(1).bound(nb1,nb2).asp_latu];
            blk(1).aline_latd = [blk(1).aline_latd; blk(1).bound(nb1,nb2).asp_latd];
          end
          % Indexing trimesh with each strip
          hx_d = (blk(1).bound(nb1,nb2).asp_xd_interp(1:end-1) + blk(1).bound(nb1,nb2).asp_xd_interp(2:end)) ./ 2;
          hx_u = (blk(1).bound(nb1,nb2).asp_xu_interp(1:end-1) + blk(1).bound(nb1,nb2).asp_xu_interp(2:end)) ./ 2;
          hy_d = (blk(1).bound(nb1,nb2).asp_yd_interp(1:end-1) + blk(1).bound(nb1,nb2).asp_yd_interp(2:end)) ./ 2;
          hy_u = (blk(1).bound(nb1,nb2).asp_yu_interp(1:end-1) + blk(1).bound(nb1,nb2).asp_yu_interp(2:end)) ./ 2;
          hx_d = [blk(1).bound(nb1,nb2).asp_xd_interp(1); hx_d; blk(1).bound(nb1,nb2).asp_xd_interp(end)];
          hx_u = [blk(1).bound(nb1,nb2).asp_xu_interp(1); hx_u; blk(1).bound(nb1,nb2).asp_xu_interp(end)];
          hy_d = [blk(1).bound(nb1,nb2).asp_yd_interp(1); hy_d; blk(1).bound(nb1,nb2).asp_yd_interp(end)];
          hy_u = [blk(1).bound(nb1,nb2).asp_yu_interp(1); hy_u; blk(1).bound(nb1,nb2).asp_yu_interp(end)];
          nstrip = size(hx_d,1)-1;
          blk(1).bound(nb1,nb2).stripid = false(nf,nstrip);
          for n = 1:nstrip
            stripx = [hx_d(n),hx_d(n+1),hx_u(n+1),hx_u(n)];
            stripy = [hy_d(n),hy_d(n+1),hy_u(n+1),hy_u(n)];
            instrip = inpolygon(tricx(mm1:mm1+nf-1),tricy(mm1:mm1+nf-1),stripx,stripy);
            blk(1).bound(nb1,nb2).stripid(:,n) = instrip;
          end
          nl = 1;
          blk(1).bound(nb1,nb2).interpid = zeros(nint*(nasp-1)+1, nasp);
          for n = 1:nasp
            if n ~= nasp
              blk(1).bound(nb1,nb2).interpid(nl:nl+nint-1,n  ) = 1-(0:nint-1)./nint;
              blk(1).bound(nb1,nb2).interpid(nl:nl+nint-1,n+1) =   (0:nint-1)./nint;
            else
              blk(1).bound(nb1,nb2).interpid(end,end) = 1;
            end
            nl = nl + nint;
          end
          blk(1).bound(nb1,nb2).naspline  = nasp;
          blk(1).bound(nb1,nb2).interpint = nint;
        else
          error(['Not found ',rwlfile]);
        end
      end
      mm1 = mm1 +   nf;
      mm3 = mm3 + 3*nf;
    end
  end
end

fprintf('=== Read Random Walk Lines === \n');
end

%% Make figures
function MakeFigs(blk,cal,bslip,obs)
% 
figure(100); clf(100)
mm3 = 1;
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    nf = size(blk(1).bound(nb1,nb2).blon,1);
    if blk(1).bound(nb1,nb2).flag2 == 1
      % Back-slip rate
      patch(blk(1).bound(nb1,nb2).blon',...
            blk(1).bound(nb1,nb2).blat',...
      sqrt(bslip(mm3:mm3+nf-1).^2+bslip(mm3+nf:mm3+2*nf-1).^2))
      % Initial patches
      for np = 1:size(blk(1).bound(nb1,nb2).patch,2)
        hold on
        plot(blk(1).bound(nb1,nb2).patch(np).lon,...
             blk(1).bound(nb1,nb2).patch(np).lat,...
             'LineWidth',3,'Color','b')
      end
      mm3 = mm3 + 3*nf;
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

%% Show results for makeing FIGURES
function MakeFigures(cha,blk,obs,rt,lo_mc,up_mc,vec,bslip,bslipl,mimean)
% Color palette(Polar)
red  = [           0: 1/32:1     ones(1,32)]';
green= [           0: 1/32:1 1-1/32:-1/32:0]';
blue = [ones(1,32) 1:-1/32:0               ]';
rwb = [red green blue];
rw  =    rwb(33:end,:);
if lo_mc==-1
  cmap=rwb;
else
  cmap=rw;
end
% Quiver plot scale factor
vfactor = 1e-2; % vector
cfactor = 5e+8; % strain

fig100 = figure(100); clf(100)
fig100.Position = [60,60,1200,900];
%---------Show kinematic coupling ------------------------
figure(100); subplot(2,2,1) % Bug to wait zero
mm1 = 1;
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    if blk(1).bound(nb1,nb2).flag2 == 1
      continue
    else
      nf = size(blk(1).bound(nb1,nb2).blon,1);
      if nf~=0
        patch(blk(1).bound(nb1,nb2).blon',blk(1).bound(nb1,nb2).blat', blk(1).bound(nb1,nb2).bdep',mean(cha.mc(mm1:mm1+nf-1,:,1),2));
        mm1 = mm1 + nf;
        hold on
      end
    end
  end
end
ax1 = gca;
ax1.CLim = [lo_mc up_mc];
colormap(ax1,cmap)
colorbar

%---------Show standard deviation for subfaults----------
figure(100); subplot(2,2,2) % Bug to wait zero
% bug to wait zero
mm1 = 1;
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    if blk(1).bound(nb1,nb2).flag2 == 1
      continue
    else
      nf = size(blk(1).bound(nb1,nb2).blon,1);
      if nf ~= 0
        patch(blk(1).bound(nb1,nb2).blon',blk(1).bound(nb1,nb2).blat',blk(1).bound(nb1,nb2).bdep',std(cha.mc(mm1:mm1+nf-1,:,1),0,2));
        mm1 = mm1 + nf;
        hold on
      end
    end
  end
end
ax2 = gca;
colormap(ax2,'parula')
colorbar

%---------Show mechanical backslip rate ------------------
mm3 = 1;
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    if blk(1).bound(nb1,nb2).flag2 == 1
      nf = size(blk(1).bound(nb1,nb2).blon,1);
      % Backslip rate (locked)
      figure(100); subplot(2,2,3); patch(blk(1).bound(nb1,nb2).blon',...
                            blk(1).bound(nb1,nb2).blat',...
                            blk(1).bound(nb1,nb2).bdep',...
                            single(sqrt(bslipl(mm3:mm3+nf-1,1).^2+bslipl(mm3+nf:mm3+2*nf-1,1).^2)~=0)); hold on
      % Backslip rate (all)
      figure(100); subplot(2,2,4); patch(blk(1).bound(nb1,nb2).blon',...
                            blk(1).bound(nb1,nb2).blat',...
                            blk(1).bound(nb1,nb2).bdep',...
                                   sqrt(bslip( mm3:mm3+nf-1,1).^2+bslip( mm3+nf:mm3+2*nf-1,1).^2)   ); hold on
      mm3 = mm3 + 3*nf;
    else
      continue
    end
  end
end
subplot(2,2,3);
ax3 = gca;
c = flipud(gray); colormap(ax3,c)
subplot(2,2,4);
ax4 = gca;
c = flipud(hot);  colormap(ax4,c);
colorbar
caxis([0, max(gather(bslipl(:,1)))*1.1]);

%---------Show 2-d histogram of sampled euler pole-------------
figure(130); clf(130)
for nb = 1:blk(1).nblock
  plot(blk(nb).lon,blk(nb).lat,'red')
  hold on
  text(mean(blk(nb).lon),mean(blk(nb).lat),int2str(nb))
  hold on
  [latp,lonp,~] = xyzp2lla(cha.mp(3.*nb-2,:,1),cha.mp(3.*nb-1,:,1),cha.mp(3.*nb,:,1));
  minlon = min(lonp); maxlon = max(lonp); 
  minlat = min(latp); maxlat = max(latp); 
  if maxlon-minlon < 0.5; binlon=[minlon maxlon]; else binlon=minlon:0.5:maxlon; end  
  if maxlat-minlat < 0.5; binlat=[minlat maxlat]; else binlat=minlat:0.5:maxlat; end  
  histogram2(lonp,latp,binlon,binlat,'normalization','probability','facecolor','flat')
  hold on
  text(double(mean(lonp)),double(mean(latp)),int2str(nb))
  hold on
end
colorbar
hold on

f140 = figure(140); clf(140)
f140.Position = [60,60,1200,900];
%---------Show obs and cal vector at sites -------------
% Color of arrows
% green : Observed velocities
% blue  : Calculated velocities
figure(140); subplot(2,2,1)
quiver(obs(1).alon,obs(1).alat,vfactor.*obs(1).evec,        vfactor.*obs(1).nvec,        'AutoScale','off','Color','green')
hold on
quiver(obs(1).alon,obs(1).alat,vfactor.*vec.sum(1:3:end,1)',vfactor.*vec.sum(2:3:end,1)','AutoScale','off','Color', 'blue')
hold on
axis([obs(1).lonmin-1,obs(1).lonmax+1,obs(1).latmin-1,obs(1).latmax+1]);
title(['Horizontal obs and cal motion (iteration number: ',num2str(rt),')']);

figure(140); subplot(2,2,2)
quiver(obs(1).alon,obs(1).alat,zeros(size(obs(1).hvec)),vfactor.*obs(1).hvec,        'AutoScale','off','Color','green')
hold on
quiver(obs(1).alon,obs(1).alat,zeros(size(obs(1).hvec)),vfactor.*vec.sum(3:3:end,1)','AutoScale','off','Color', 'blue')
hold on
axis([obs(1).lonmin-1,obs(1).lonmax+1,obs(1).latmin-1,obs(1).latmax+1]);
title(['Vertical obs and cal motion (iteration number: ',num2str(rt),')']);

%-------------------- Show rigid, kinematic, mechanical vectors------------
% Color of arrows
% black   : Rigid rotation
% red     : Elastic deformation due to mechanical coupling
% blue    : Elastic deformation due to kinematic coupling
figure(140); subplot(2,2,3)
quiver(obs(1).alon,obs(1).alat, vfactor.*vec.rig(1:3:end,1)',vfactor.*vec.rig(2:3:end,1)','AutoScale','off','Color','k')
hold on
quiver(obs(1).alon,obs(1).alat, vfactor.*vec.kin(1:3:end,1)',vfactor.*vec.kin(2:3:end,1)','AutoScale','off','Color','b')
hold on
quiver(obs(1).alon,obs(1).alat, vfactor.*vec.mec(1:3:end,1)',vfactor.*vec.mec(2:3:end,1)','AutoScale','off','Color','r')
axis([obs(1).lonmin-1,obs(1).lonmax+1,obs(1).latmin-1,obs(1).latmax+1]);
title(['Horizontal rigid and elastic motion (iteration number: ',num2str(rt),')']);

figure(140); subplot(2,2,4)
quiver(obs(1).alon,obs(1).alat, zeros(size(vec.rig(3:3:end,1)))',vfactor.*vec.rig(3:3:end,1)','AutoScale','off','Color','k')
hold on
quiver(obs(1).alon,obs(1).alat, zeros(size(vec.kin(3:3:end,1)))',vfactor.*vec.kin(3:3:end,1)','AutoScale','off','Color','b')
hold on
quiver(obs(1).alon,obs(1).alat, zeros(size(vec.mec(3:3:end,1)))',vfactor.*vec.mec(3:3:end,1)','AutoScale','off','Color','r')
axis([obs(1).lonmin-1,obs(1).lonmax+1,obs(1).latmin-1,obs(1).latmax+1]);
title(['Vertical rigid and elastic motion (iteration number: ',num2str(rt),')']);

%------------------------Show principal strain-----------------------------
% Color of arrows
% cyan    : Principal strain 1, maximum
% magenta : Principal strain 2, minimum
figure(160); clf(160)
for nb = 1:blk(1).nblock
  hold on; plot(blk(nb).lon,blk(nb).lat,'red')
  hold on; text(mean(blk(nb).lon),mean(blk(nb).lat),int2str(nb),'color','r')
  e=[mimean(3*nb-2,1) mimean(3*nb-1,1);...
     mimean(3*nb-1,1) mimean(3*nb  ,1)];
  [eigv,eigd] = eig(e);
  e1 = eigd(1,1); e2 = eigd(2,2);
  v1 = eigv(:,1); v2 = eigv(:,2);
  if e1>=0; c1='c'; else; c1='m'; end
  if e2>=0; c2='c'; else; c2='m'; end
  figure(160)
  hold on; quiver(blk(nb).loninter,blk(nb).latinter,cfactor* e1*v1(1),cfactor* e1*v1(2),'AutoScale','off','Color',c1,'showarrowhead','off','linewidth',1);
  hold on; quiver(blk(nb).loninter,blk(nb).latinter,cfactor*-e1*v1(1),cfactor*-e1*v1(2),'AutoScale','off','Color',c1,'showarrowhead','off','linewidth',1);
  hold on; quiver(blk(nb).loninter,blk(nb).latinter,cfactor* e2*v2(1),cfactor* e2*v2(2),'AutoScale','off','Color',c2,'showarrowhead','off','linewidth',1);
  hold on; quiver(blk(nb).loninter,blk(nb).latinter,cfactor*-e2*v2(1),cfactor*-e2*v2(2),'AutoScale','off','Color',c2,'showarrowhead','off','linewidth',1);
  hold on; plot(blk(nb).loninter,blk(nb).latinter,'.k','markersize',5)
  hold on; text(blk(nb).loninter,blk(nb).latinter,int2str(nb),'color','k')
end
axis([obs(1).lonmin-1,obs(1).lonmax+1,obs(1).latmin-1,obs(1).latmax+1]);
title(['principle strain (iteration number: ',num2str(rt),')']);

% debug----------
drawnow
end

%% Wait GPU program
function free=waitGPU(varargin)
% Coded by Zhang Xuelei
a=true;
d=gpuDevice;
if isempty(varargin)
  limit=30;
else
  limit=varargin{1};
end
% tic
while a
  if limit<=100
    free=d.FreeMemory/d.TotalMemory*100;
    if free>limit; break; end
    pause(0.5);
  elseif limit>100
    free=d.FreeMemory;
    if free>limit; break; end
    pause(0.5)
  end
end
% waittime=toc;
% if waittime>0.5
%     disp(['waiting time' num2str(waittime)])
% end
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

%% Heaviside step function
function y = Heaviside(x)
% Calculate Heveaside step function
%   y = Heaviside(x)
%   y and x are possible to be either scaler or vector.
y = 1 .* (x >= 0);
end

%% Interp linearly for some points
function y = linspace2(x,nint)
% Calculate linspace for vector data
% y = linspace2(x, nint) x and y, vector; nint, interp interval with
tmpx = x;
if size(tmpx,1) > size(tmpx,2)
  x = x';
end
x1 = x(1:end-1);
x2 = x(2:end  );
xdiff = x2 - x1;
n1 = nint-1;
y = x1 + (xdiff./nint) .* (0:n1)';
y = [y(:); x(end)];
if size(tmpx,1) < size(tmpx,2)
  y = y';
end
end