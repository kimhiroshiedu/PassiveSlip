function CalcSlipArea(savefolder,blk,obs,G,d,tcha,lock,T,event,threratio)
% Load blk, grn, obs, tcha and lock.mat files before running this script.
% event
% segment_number1  max_slip1
% segment_number2  max_slip2
%  :                  :
%  :                  :
G(1).tb_kin = full(G(1).tb_kin);
G(1).tb_mec = full(G(1).tb_mec);

idl   = zeros(  size(tcha.aveaid,1),1);
slipl = zeros(3*size(tcha.aveaid,1),1);
for eq = event'
  [s] = calc_slip_asperity(tcha,G,d,lock,T,eq,-1);
  [rt,lock] = est_reccurence_time(blk,lock,s,eq);
  [s] = calc_slip_asperity(tcha,G,d,lock,T,eq,rt);
  idl = idl | lock(eq(1)).idl50;
  slipl = slipl + s.slipl50;
end
[s]   = calc_slip_total(d,G,slipl,idl);
[s]   = calc_seismic_moment(blk,lock,event,threratio,s);
[s,v] = calc_surface_velocity(blk,obs,G,d,s);
save_result(obs,blk,tcha,lock,s,v,T,event,savefolder)
end

function [s] = calc_slip_asperity(tcha,G,d,lock,T,eq,rt)
event = eq(1);
sslip = eq(2);
mpmean = tcha.avepol(:,:,T);
idl50 = logical(d(1).maid *  lock(event).idl50);

% Calculate back-slip on locked and creeping patches.
rel_mec    = -rt .* (G(1).tb_mec * mpmean) .* d(1).cfinv_mec;
slipl50    = rel_mec .* idl50;

s.idl50    = idl50;
s.rel_mec  = rel_mec;
s.slipl50  = slipl50;
end

function [rt,lock] = est_reccurence_time(blk,lock,s,eq)
event = eq(1);
sslip = eq(2);
idl = lock(event).idl50;
lock(event).rt = [];
mm3m = 1;
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    nf = size(blk(1).bound(nb1,nb2).blon,1);
    if nf ~= 0
      if blk(1).bound(nb1,nb2).flag2 == 1
        slip50_st = s.slipl50(mm3m     :mm3m+  nf-1);
        slip50_dp = s.slipl50(mm3m+  nf:mm3m+2*nf-1);
        slip50_ts = s.slipl50(mm3m+2*nf:mm3m+3*nf-1);
        bslip50_s  = sqrt(slip50_st.^2 + slip50_dp.^2 + slip50_ts.^2);
        rt = sslip ./ bslip50_s;
        lock(event).rt = [lock(event).rt; rt];
        mm3m = mm3m + 3*nf;
      end
    end
  end
end
rt = max(lock(event).rt(idl));
end

function [s] = calc_slip_total(d,G,slipl,idl)
idl50 = logical(d(1).maid *  idl);
idc50 = logical(d(1).maid * ~idl);
slip50 = slipl;
slip50(idc50) = -G(1).s(idc50,idc50) \ (G(1).s(idc50,idl50) * slipl(idl50));
s.idl   = idl50;
s.idc   = idc50;
s.slipl = slipl;
s.slip  = slip50;
end

function [s] = calc_seismic_moment(blk,lock,event,threratio,s)
mu = 40; % (GPa)
uf = 1e12; % unit translation factor

combibo = [];
for eq = event'
  combibo = [combibo; lock(eq(1)).boundary];
end
combibo = unique(combibo,'rows');

mm3m = 1;
maxslipl = 0;
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    nf = size(blk(1).bound(nb1,nb2).blon,1);
    if nf ~= 0
      if blk(1).bound(nb1,nb2).flag2 == 1
        slipl = sqrt(s.slipl(mm3m     :mm3m+  nf-1).^2 ...
                   + s.slipl(mm3m+  nf:mm3m+2*nf-1).^2 ...
                   + s.slipl(mm3m+2*nf:mm3m+3*nf-1).^2);
        maxslipl = max(maxslipl,max(slipl));
        tmp(1).bound(nb1,nb2).slip_st  = s.slip(mm3m     :mm3m+  nf-1);
        tmp(1).bound(nb1,nb2).slip_dp  = s.slip(mm3m+  nf:mm3m+2*nf-1);
        tmp(1).bound(nb1,nb2).slip_ts  = s.slip(mm3m+2*nf:mm3m+3*nf-1);
        tmp(1).bound(nb1,nb2).triarea  = lock(1).bound(nb1,nb2).area;
        mm3m = mm3m + 3*nf;
      end
    end
  end
end
Mo     =  0;
Mo_cut =  0;
S = 0;
D = 0;
cutoff_total = 0;
for bo = combibo'
  nb1 = min(bo);
  nb2 = max(bo);
  slip_st = tmp(1).bound(nb1,nb2).slip_st;
  slip_dp = tmp(1).bound(nb1,nb2).slip_dp;
  slip_ts = tmp(1).bound(nb1,nb2).slip_ts;
  area    = tmp(1).bound(nb1,nb2).triarea;
  slip    = sqrt(slip_st.^2 + slip_dp.^2 + slip_ts.^2);
  cutoff  = slip >= threratio*maxslipl;
  sliparea= cutoff .* area;
  Mo_i     = uf .* mu .*            slip'  * area;
  Mo_i_cut = uf .* mu .* (cutoff .* slip)' * area;
  S = S + sum(sliparea);
  D = D + cutoff' * slip;
  cutoff_total = cutoff_total + sum(cutoff);
  Mo     = Mo + Mo_i;
  Mo_cut = Mo_cut + Mo_i_cut;
  s(1).bound(nb1,nb2).slip_st = tmp(1).bound(nb1,nb2).slip_st;
  s(1).bound(nb1,nb2).slip_dp = tmp(1).bound(nb1,nb2).slip_dp;
  s(1).bound(nb1,nb2).slip_ts = tmp(1).bound(nb1,nb2).slip_ts;
  s(1).bound(nb1,nb2).slip    = slip;
  s(1).bound(nb1,nb2).sliparea= sliparea;
  s(1).bound(nb1,nb2).cutoff  = cutoff;
  s(1).bound(nb1,nb2).Mo     = Mo_i;
  s(1).bound(nb1,nb2).Mo_cut = Mo_i_cut;
  s(1).bound(nb1,nb2).Dmean = mean(slip(cutoff));
  s(1).bound(nb1,nb2).S     = sum(sliparea);
end

s(1).mu = mu;
s(1).threratio = threratio;
s(1).threslip  = threratio * maxslipl;
s(1).Mo = Mo;
s(1).Mo_cut = Mo_cut;
s(1).S = S;
s(1).D_mean = D / cutoff_total;
end

function [s,v] = calc_surface_velocity(blk,obs,G,d,s)
s(1).idcutoff = false(size(s.slip,1)/3,1);
mm1m = 1;
mm3m = 1;
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    nf = size(blk(1).bound(nb1,nb2).blon,1);
    if nf ~= 0
      if blk(1).bound(nb1,nb2).flag2 == 1
        slip = sqrt(s.slip(mm3m     :mm3m+  nf-1).^2 ...
                  + s.slip(mm3m+  nf:mm3m+2*nf-1).^2 ...
                  + s.slip(mm3m+2*nf:mm3m+3*nf-1).^2);
        s(1).idcutoff(mm1m:mm1m+nf-1) = slip >= s(1).threslip;
        mm1m = mm1m +   nf;
        mm3m = mm3m + 3*nf;
      end
    end
  end
end
v.coseis        = G(1).c_mec *  s.slip;
v.coseis_cutoff = G(1).c_mec * (s.slip .* (d(1).maid * s(1).idcutoff));
v.obs = reshape([obs(1).evec; obs(1).nvec; obs(1).hvec],3*obs(1).nobs,1);
end

function save_result(obs,blk,tcha,lock,s,v,T,event,savefolder)
if exist(savefolder,'dir')~=7; mkdir(savefolder); end
% Save coseismic slip magnitude
mm1  = 1;
mm3m = 1;
mm1m = 1;
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    nf = size(blk(1).bound(nb1,nb2).blon,1);
    if nf ~= 0
      fltnum = mm1:mm1+nf-1;
      if blk(1).bound(nb1,nb2).flag2 == 1
        clon = mean(blk(1).bound(nb1,nb2).blon,2);
        clat = mean(blk(1).bound(nb1,nb2).blat,2);
        cdep = mean(blk(1).bound(nb1,nb2).bdep,2);
        slip_st  = s.slip(mm3m     :mm3m+  nf-1);
        slip_dp  = s.slip(mm3m+  nf:mm3m+2*nf-1);
        slip_ts  = s.slip(mm3m+2*nf:mm3m+3*nf-1);
        id_lock  = s.idl( mm3m     :mm3m+  nf-1);
        id_slip  = sqrt(slip_st.^2+slip_dp.^2+slip_ts.^2) >= s(1).threslip;
        pc = tcha.aveaid(mm1m:mm1m+nf-1,:,T);
        outdata = [fltnum',...
            blk(1).bound(nb1,nb2).blon,...
            blk(1).bound(nb1,nb2).blat,...
            blk(1).bound(nb1,nb2).bdep,...
            clon, clat, cdep, pc,id_lock,id_slip,...
            slip_st.*1e-3, slip_dp.*1e-3, slip_ts.*1e-3];
        file = fullfile(savefolder,['trimec_',num2str(nb1),'_',num2str(nb2),'.txt']);
        fid  = fopen(file,'wt');
        fprintf(fid,'#    1        2        3        4       5       6       7       8       9      10       11      12      13   14   15     16      17      18      19\n');
        fprintf(fid,'#  tri     lon1     lon2     lon3    lat1    lat2    lat3    dep1    dep2    dep3     clon    clat    cdep  P_l id_l slipid slip_st slip_dp slip_ts [m]\n');
        fprintf(fid,'%6i %8.3f %8.3f %8.3f %7.3f %7.3f %7.3f %7.2f %7.2f %7.2f %8.3f %7.3f %7.2f %4.2f %4i %6i %7.3f %7.3f %7.3f \n',outdata');
        fclose(fid);
        mm1m = mm1m +   nf;
        mm3m = mm3m + 3*nf;        
      end
      mm1 = mm1 + nf;
    end
  end
end
s.event = event;
save(fullfile(savefolder,'slip'),'s','-v7.3');

% Save released moment
combibo = [];
for eq = event'
  combibo = [combibo; lock(eq(1)).boundary];
end
combibo = unique(combibo,'rows');
fid = fopen(fullfile(savefolder,'released_moment.txt'),'wt');
fprintf(fid,'# b1 b2 meanD(m)    S(km^2)     Mo(Nm) Mo_thr(Nm)\n');
for bo = combibo'
  nb1 = min(bo);
  nb2 = max(bo);
  fprintf(fid,' %2i %2i %8.2e %10.2e %10.2e %10.2e\n',...
      bo(1),bo(2),...
      s(1).bound(nb1,nb2).Dmean*1e-3,...
      s(1).bound(nb1,nb2).S,...
      s(1).bound(nb1,nb2).Mo,...
      s(1).bound(nb1,nb2).Mo_cut);
end
fprintf(fid,'# Total moment\n');
fprintf(fid,'       %8.2e %10.2e %10.2e %10.2e\n',...
    s(1).D_mean*1e-3,s(1).S,s(1).Mo,s(1).Mo_cut);
fprintf(fid,'# Ruccurence time (yr)\n');
for eq = event'
  fprintf(fid,'%15s %5.0f\n',lock(eq(1)).name,max(lock(eq(1)).rt(lock(eq(1)).idl50)));
end
fclose(fid);

% Save vectors
file = fullfile(savefolder,'cal_vector.txt');
ve = v.coseis(1:3:end);
vn = v.coseis(2:3:end);
vu = v.coseis(3:3:end);
ve_cut = v.coseis_cutoff(1:3:end);
vn_cut = v.coseis_cutoff(2:3:end);
vu_cut = v.coseis_cutoff(3:3:end);
fid = fopen(file,'wt');
fprintf(fid,'#   site      lon     lat        hei     ve     vn     vu  vecut  vncut  vucut (mm/yr)\n');
for nob = 1:obs(1).nobs
  fprintf(fid,'%8s %8.3f %7.3f %10.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n',...
      obs(1).name{nob},obs(1).alon(nob),obs(1).alat(nob),obs(1).ahig(nob),...
      ve(nob).*1e-3,vn(nob).*1e-3,vu(nob).*1e-3,...
      ve_cut(nob).*1e-3,vn_cut(nob).*1e-3,vu_cut(nob).*1e-3);
end
fclose(fid);

end
