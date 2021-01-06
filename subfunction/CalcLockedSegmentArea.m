function CalcLockedSegmentArea(savefolder,prm,blk,obs,G,d,tcha,patchfolder,T)
% T, Temperature of replica
% rt, reccurence time to calculate slip magnitude (positive, coseismic; negative, interseismic)
G(1).tb_kin = full(G(1).tb_kin);
G(1).tb_mec = full(G(1).tb_mec);
[blk,lock] = ReadLockedPatch(blk,patchfolder);
[s,~] = CalcOptimumValue(prm,obs,tcha,G,d,T);
SaveAsperitySegmentArea(savefolder,blk,obs,tcha,lock,T,s);
end

%%
function SaveAsperitySegmentArea(savefolder,blk,obs,tcha,lock,T,s)
if exist(savefolder,'dir')~=7; mkdir(savefolder); end
alat0 = mean(obs(1).alat,2);
alon0 = mean(obs(1).alon,2);
mm1m = 1;
mm3m = 1;
idl50 = zeros(size(tcha.aveaid,1),1);
tricx = zeros(size(tcha.aveaid,1),1);
tricy = zeros(size(tcha.aveaid,1),1);
triarea = zeros(size(tcha.aveaid,1),1);
triline = zeros(size(tcha.aveaid,1),1);
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    nf = size(blk(1).bound(nb1,nb2).blon,1);
    if nf ~= 0
      if blk(1).bound(nb1,nb2).flag2 == 1
        [trix,triy] = PLTXY(blk(1).bound(nb1,nb2).blat,blk(1).bound(nb1,nb2).blon,alat0,alon0);
        triz = blk(1).bound(nb1,nb2).bdep;
        [cx,cy] = PLTXY(mean(blk(1).bound(nb1,nb2).blat,2),mean(blk(1).bound(nb1,nb2).blon,2),alat0,alon0);
        area = zeros(nf,1);
        for ntri = 1:size(blk(1).bound(nb1,nb2).blat,1)
          area(ntri) = triangle_area([trix(ntri,:)', triy(ntri,:)', triz(ntri,:)']);
        end
        blk(1).bound(nb1,nb2).triid = false(size(tcha.aveaid,1),1);
        blk(1).bound(nb1,nb2).triid(mm1m:mm1m+nf-1) = true;
        tricx(mm1m:mm1m+nf-1) = cx;
        tricy(mm1m:mm1m+nf-1) = cy;
        triarea(mm1m:mm1m+nf-1) = area;
        triline(mm1m:mm1m+nf-1) = sqrt((4/sqrt(3)).*area);
        lock(1).bound(nb1,nb2).area = triarea(mm1m:mm1m+nf-1);
        lock(1).bound(nb1,nb2).line = triline(mm1m:mm1m+nf-1);
        idl50(mm1m:mm1m+nf-1) = s.idl50(mm3m:mm3m+nf-1);
        mm1m = mm1m +   nf;
        mm3m = mm3m + 3*nf;
      end
    end
  end
end

fid = fopen(fullfile(savefolder,'locksegments_area.txt'),'wt');
fprintf(fid,'#     segment   area(km^2)\n');
for nlock = 1:size(lock,2)
  id_lock = zeros(size(tcha.aveaid,1),1);
  [edgex,edgey] = PLTXY(lock(nlock).lat,lock(nlock).lon,alat0,alon0);
  id_segment = inpolygon(tricx,tricy,edgex,edgey);
  for nbound = 1:size(lock(nlock).boundary,1)
    nb1 = min(lock(nlock).boundary(nbound,:));
    nb2 = max(lock(nlock).boundary(nbound,:));
    id = blk(1).bound(nb1,nb2).triid .* id_segment .* idl50;
    id_lock = id_lock | id;
  end
  lock(nlock).area = triarea' * id_lock;
  lock(nlock).idl50 = id_lock;
  fprintf(fid,'%15s %10.2f\n',lock(nlock).name,lock(nlock).area);
end
save(fullfile(savefolder,'lock'),'lock');
fprintf('=== Pass SaveAsperitySegmentArea === \n');
end

%% Save vector of seismic cycle
function SaveVectors(folder,obs,vec,T)
savedir = fullfile(folder,'vector');
if exist(savedir) ~=7; mkdir(savedir); end
% Save vmec
file = fullfile(savedir,'mec_vector.txt');
savevector(file,obs,vec.mec,vec.mec50);
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

%% Calc optimum value from tcha.mat
function [s,vec] = CalcOptimumValue(prm,obs,tcha,G,d,T)
mpmean = tcha.avepol(:,:,T);
mcmean = tcha.aveflt(:,:,T);
mamean = tcha.aveasp(:,:,T);
mimean = tcha.aveine(:,:,T);
idmean = tcha.aveaid(:,:,T);

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
% load('Result_xps/Test_06/locksegments/kimura2019/highsdr_KH2019.mat')  % High-SDR (Kimura et al., 2019)
load('Result_red/Test_27/locksegments/hashimoto2012/highsdr_HC2012.mat')  % High-SDR (Hashimoto et al., 2012)
% load('lock_id_1968eq.mat')  % 1968 asperity
idmean(c_idin) = 1;
idmean(~c_idin) = 0;
% -----
idl25 = logical(d(1).maid *  (idmean>=0.25));
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

fprintf('=== Pass CalcOptimumValue === \n');

end

%% Read locked patches
function [blk,lock] = ReadLockedPatch(blk,patchfolder)
%    Test version coded by H. Kimura 2019/1/29
% Revised version coded by H. Kimura 2019/2/5

patchfile = fullfile(patchfolder,'locksegments.txt');
fid       = fopen(patchfile,'r');
if fid > 0
  np = 0;
  while 1
    n  = 0;
    tline = fgetl(fid);
    if ~ischar(tline)
      break
    elseif tline(1) == '#'
      continue
    elseif tline(1) == '@'
      lstr = strsplit(tline);
      share = zeros((size(lstr,2)-1)/2,2);
      for nshare = 1:size(share,1)
        share(nshare,:) = [str2num(lstr{1+2*nshare-1}),str2num(lstr{1+2*nshare})];
      end
      while 1
        tline = fgetl(fid);
        if ~ischar(tline) || tline(1) == '#'; break; end
        if tline(1) == '>'
          np = np+1;
          lstr = strsplit(tline);
          if size(lstr,2) > 1
            lock(np).name = strjoin(lstr(2:end));
          else
            lock(np).name = num2str(np);
          end
          lock(np).boundary = share;
          n  = 0;
          continue
        elseif tline(1) ~= '>' && tline(1) ~= '@'
          n   = n+1;
          tmp = strsplit(strtrim(tline));
          lock(np).lon(n) = str2double(cellstr(tmp(1)));
          lock(np).lat(n) = str2double(cellstr(tmp(2)));
        elseif tline(1) == '@'
          error('Leave one line between segment combinations')
        end
      end
    end
  end
end
fprintf('=== Pass ReadLockedPatches === \n');
end

%% Subroutines
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

function y = Heaviside(x)
% Calculate Heveaside step function
%   y = Heaviside(x)
%   y and x are possible to be either scaler or vector.
y = 1 .* (x >= 0);
end
