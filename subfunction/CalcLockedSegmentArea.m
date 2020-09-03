function CalcLockedSegmentArea(savefolder,prm,blk,obs,G,d,tcha,patchfolder,T)
% T, Temperature of replica
% rt, reccurence time to calculate slip magnitude (positive, coseismic; negative, interseismic)
G(1).tb_kin = full(G(1).tb_kin);
G(1).tb_mec = full(G(1).tb_mec);
blk = ReadLockedPatch(blk,patchfolder);
[s,~] = CalcOptimumValue(prm,obs,tcha,G,d,T);
SaveAsperitySegmentArea(savefolder,blk,obs,tcha,T,s);
end

%%
function SaveAsperitySegmentArea(folder,blk,obs,tcha,T,s)
savefolder = fullfile(folder,['T_',num2str(T,'%02i')],'locksegments');
if exist(savefolder,'dir')~=7; mkdir(savefolder); end
alat0 = mean(obs(1).alat,2);
alon0 = mean(obs(1).alon,2);
mm1m = 1;
mm3m = 1;
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
      triasp(1).bound(nb1,nb2).area = area;
      triasp(1).bound(nb1,nb2).line = sqrt((4/sqrt(3)).*area);
      if blk(1).bound(nb1,nb2).flag2 == 1
        if blk(1).bound(nb1,nb2).segmentid == 1
          fid = fopen(fullfile(savefolder,['lockingarea_',num2str(nb1,'%i'),'_',num2str(nb2,'%i'),'.txt']),'wt');
          fprintf(fid,'#     segment   area(km^2)\n');
          for nseg = 1:size(blk(1).bound(nb1,nb2).segment,2)
            triasp(1).bound(nb1,nb2).segment(nseg).triid = false(size(tcha.aveaid,1),1);
            [edgex,edgey] = PLTXY(blk(1).bound(nb1,nb2).segment(nseg).lat,blk(1).bound(nb1,nb2).segment(nseg).lon,alat0,alon0);
            segid = inpolygon(tricx,tricy,edgex,edgey);
            triasp(1).bound(nb1,nb2).segment(nseg).name = blk(1).bound(nb1,nb2).segment(nseg).name;
            triasp(1).bound(nb1,nb2).segment(nseg).triid(mm1m:mm1m+nf-1) = segid .* s.idl50(mm3m:mm3m+nf-1);
            triasp(1).bound(nb1,nb2).segment(nseg).lockarea = segid'.*area' * s.idl50(mm3m:mm3m+nf-1,:);
            fprintf(fid,'%15s %10.2f\n',...
                triasp(1).bound(nb1,nb2).segment(nseg).name,triasp(1).bound(nb1,nb2).segment(nseg).lockarea);
          end
          fclose(fid);
        else
          triasp(1).bound(nb1,nb2).smparea = blk(1).bound(nb1,nb2).triarea' * tcha.smpaid(mm1m:mm1m+nf-1,:);
        end
        mm1m = mm1m +   nf;
        mm3m = mm3m + 3*nf;  
      end
    end
  end
end
save(fullfile(savefolder,'triasp'),'triasp');
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
bslip25(idc25) = -G(1).s(idc25,idc25) \ (G(1).s(idc25,idl25) * bslip(idl25));
bslip50(idc50) = -G(1).s(idc50,idc50) \ (G(1).s(idc50,idl50) * bslip(idl50));
bslip75(idc75) = -G(1).s(idc75,idc75) \ (G(1).s(idc75,idl75) * bslip(idl75));
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

%% Read locked patches
function [blk] = ReadLockedPatch(blk,patchfolder)
%    Test version coded by H. Kimura 2019/1/29
% Revised version coded by H. Kimura 2019/2/5

for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    if blk(1).bound(nb1,nb2).flag2 == 1
      blk(1).bound(nb1,nb2).segmentid = 0;
      patchfile = fullfile(patchfolder,['locksegments_',num2str(nb1),'_',num2str(nb2),'.txt']);
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
            blk(1).bound(nb1,nb2).segment(np).lon(n) = str2double(cellstr(tmp(1)));
            blk(1).bound(nb1,nb2).segment(np).lat(n) = str2double(cellstr(tmp(2)));
          elseif tline(1) == '>'
            np = np+1;
            lstr = strsplit(tline);
            if size(lstr,2) > 1
              blk(1).bound(nb1,nb2).segment(np).name = strjoin(lstr(2:end));
            else
              blk(1).bound(nb1,nb2).segment(np).name = num2str(np);
            end
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