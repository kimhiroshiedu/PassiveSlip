%% Main part
function out_passiveanalysis_v1(folder)
result_dir = fullfile(pwd,folder);
fprintf('Now loading %s ...',fullfile(result_dir,'/obs.mat'))
load(fullfile(result_dir,'/obs.mat'));fprintf('load\n')
fprintf('Now loading %s ...',fullfile(result_dir,'/blk.mat'))
load(fullfile(result_dir,'/blk.mat'));fprintf('load\n')
fprintf('Now loading %s ...',fullfile(result_dir,'/cal.mat'))
load(fullfile(result_dir,'/cal.mat'));fprintf('load\n')
fprintf('Now loading %s ...',fullfile(result_dir,'/grn.mat'))
load(fullfile(result_dir,'/grn.mat'));fprintf('load\n')

ExportSlipDeficit(result_dir,cal,blk,d);
ExportAsperities(folder,blk);
ExportVectors(result_dir,obs,cal);

fprintf('Files exported! \n');
end

%% Save slip deficit to txt files
function ExportSlipDeficit(folder,cal,blk,d)
mc = 1;
mr = 1;
dir_out=[folder,'/backslip'];
exid=exist(dir_out);
if exid~=7; mkdir(dir_out); end

for nb1=1:blk(1).nblock
  for nb2=nb1+1:blk(1).nblock
    nf=size(blk(1).bound(nb1,nb2).blon,1);
    if nf~=0
      fid     = fopen([dir_out,'/bs_',num2str(nb1),'_',num2str(nb2),'.txt'],'w');
      clon    = mean(blk(1).bound(nb1,nb2).blon,2);
      clat    = mean(blk(1).bound(nb1,nb2).blat,2);
      cdep    = mean(blk(1).bound(nb1,nb2).bdep,2);
      sdr     = sqrt( cal.slip(mc     :mc+  nf-1).^2 ...
                     + cal.slip(mc+  nf:mc+2*nf-1).^2 ...
                     + cal.slip(mc+2*nf:mc+3*nf-1).^2 );
      sdr_int = sdr.*d(1).id_lock(mc:mc+nf-1);
      flt_id  = mr:mr+nf-1;
      
      outdata = [flt_id' ...
          blk(1).bound(nb1,nb2).blon ...
          blk(1).bound(nb1,nb2).blat ...
          blk(1).bound(nb1,nb2).bdep ...
          clon clat cdep ...
          sdr sdr_int];
      
      fprintf(fid,'# %6s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %10s %10s \n',...
                      'tri_no','lon1','lon2','lon3','lat1','lat2','lat3','dep1','dep2','dep3',...
                      'c_lon','c_lat','c_dep','sdr[mm/yr]','sdr_int');
      fprintf(fid,'%8d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %10.4f %10.4f \n',outdata');
      fclose(fid);
      mr = mr+  nf;
      mc = mc+3*nf;
    end
  end
end

end

%% Export asperities
function ExportAsperities(folder,blk)
dir_out=[folder,'/backslip'];
exid=exist(dir_out);
if exid~=7; mkdir(dir_out); end
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    nf = size(blk(1).bound(nb1,nb2).blon,1);
    if nf ~= 0
      fid = fopen([dir_out,'/patchb_',num2str(nb1),'_',num2str(nb2),'.txt'],'w');
      for np = 1:size(blk(1).bound(nb1,nb2).patch,2)
        fprintf(fid,'> \n');
        patch = [blk(1).bound(nb1,nb2).patch(np).lon; ...
                 blk(1).bound(nb1,nb2).patch(np).lat];
        fprintf(fid,'%10.4f %10.4f \n',patch);
      end
      fclose(fid);
    end
  end
end

end

%% Export vectors of calculation, observation, residual vector
function ExportVectors(folder,obs,cal)
% 
% Calculate vector based on sampled parameter and green function
calvec.rig=cal.rig;
calvec.ela=cal.ela;
calvec.ine=cal.ine;
calvec.sum=cal.smp;
% Calculate residual vectors
obsv=[obs(1).evec;obs(1).nvec;obs(1).hvec];
calv=[cal.smp(1:3:end)';cal.smp(2:3:end)';cal.smp(3:3:end)'];
resvec.sum=obsv-calv;
% 
% Save to txt files
odir=[folder,'/vector'];
exid=exist(odir);
if exid~=7; mkdir(odir); end
obsvec=[obs(1).alon;obs(1).alat;obs(1).evec;obs(1).nvec;obs(1).hvec];
calvec.rig=[obs(1).alon;obs(1).alat;calvec.rig(1:3:end)';calvec.rig(2:3:end)';calvec.rig(3:3:end)'];
calvec.ela=[obs(1).alon;obs(1).alat;calvec.ela(1:3:end)';calvec.ela(2:3:end)';calvec.ela(3:3:end)'];
calvec.sum=[obs(1).alon;obs(1).alat;calvec.sum(1:3:end)';calvec.sum(2:3:end)';calvec.sum(3:3:end)'];
resvec.sum=[obs(1).alon;obs(1).alat;resvec.sum];
fid=fopen([odir,'/obs_vector.txt'],'w');
fprintf(fid,'%f %f %f %f %f\n',obsvec);
fclose(fid);
fid=fopen([odir,'/cal_vector.txt'],'w');
fprintf(fid,'%f %f %f %f %f\n',calvec.sum);
fclose(fid);
fid=fopen([odir,'/cal_vector_rig_site.txt'],'w');
fprintf(fid,'%f %f %f %f %f\n',calvec.rig);
fclose(fid);
fid=fopen([odir,'/cal_vector_ela_site.txt'],'w');
fprintf(fid,'%f %f %f %f %f\n',calvec.ela);
fclose(fid);
fid=fopen([odir,'/res_vector.txt'],'w');
fprintf(fid,'%f %f %f %f %f\n',resvec.sum);
fclose(fid);
% 
end
