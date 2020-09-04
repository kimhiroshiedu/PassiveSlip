function CalcSlipArea(savefolder,prm,blk,obs,G,d,tcha,patchfolder,event,T)
% event
% segment_number1  max_slip1
% segment_number2  max_slip2
%  :                  :
%  :                  :
[s,vec] = CalcOptimumValue(prm,obs,tcha,G,d,T);
SaveSlipArea(folder,blk,obs,tcha,lock,T,s);
SaveVectors(folder,obs,vec,T);
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

fprintf('=== Pass CalcOptimumValue === \n');

end

%%
function SaveSlipArea(folder,blk,obs,tcha,lock,T,s)
savefolder = fullfile(folder,['T_',num2str(T,'%02i')],'locksegments');
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

