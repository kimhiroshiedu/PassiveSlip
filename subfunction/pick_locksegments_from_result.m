% before using this script, load blk, tcha, grn,
G(1).tb_kin = full(G(1).tb_kin);
G(1).tb_mec = full(G(1).tb_mec);
%% Calculate back-slip on locked and creeping patches.
mpmean = tcha.avepol(:,:,T);
mamean = tcha.aveasp(:,:,T);
idmean = tcha.aveaid(:,:,T);

Hu = Heaviside(G(1).zulim - G(1).zc);
Hd = Heaviside(G(1).zdlim - G(1).zc);
Hlim = Hd - Hu;

idl1 = ((Heaviside(G(1).zd*mamean-G(1).zc) - Heaviside(G(1).zu*mamean-G(1).zc)) > 0) .* Hlim;
idl = logical(d(1).maid *  idl1);
idc = logical(d(1).maid * ~idl1);
idl50 = logical(d(1).maid *  (idmean>=0.50));
idc50 = logical(d(1).maid * ~(idmean>=0.50));

rel_mec    = (G(1).tb_mec * mpmean) .* d(1).cfinv_mec;
bslip      = rel_mec .* idl;
bslip50    = rel_mec .* idl50;
bslipl     = bslip;
bslipl50   = bslip50;

bslip(idc) = -G(1).s(idc,idc) \ (G(1).s(idc,idl) * bslip(idl));
bslip50(idc50) = -G(1).s(idc50,idc50) \ (G(1).s(idc50,idl50) * bslip50(idl50));

%% Make figure showing mesh, coast and block lines, and slip rate
figure(100); clf(100)
mm3 = 1;
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    if blk(1).bound(nb1,nb2).flag2 == 1
      nf = size(blk(1).bound(nb1,nb2).blon,1);
      % Backslip rate (locked)
      patch(blk(1).bound(nb1,nb2).blon',blk(1).bound(nb1,nb2).blat',...
            single(sqrt(bslipl50(mm3:mm3+nf-1,1).^2+bslipl50(mm3+nf:mm3+2*nf-1,1).^2)~=0));
      hold on
      % Backslip rate (total)
      % patch(blk(1).bound(nb1,nb2).blon',blk(1).bound(nb1,nb2).blat',...
      %       single(sqrt(bslip50(mm3:mm3+nf-1,1).^2+bslip50(mm3+nf:mm3+2*nf-1,1).^2)~=0));
      mm3 = mm3 + 3*nf;
    else
      continue
    end
  end
end
colormap(flipud(gray))

% Coast line
latlim   = [ 10  50];
lonlim   = [120 150];
filename = gunzip('gshhs_i.b.gz', tempdir);
japan    = gshhs(filename{1},latlim,lonlim);
geoshow([japan.Lat], [japan.Lon], 'displaytype','line','color','c')
ax = gca;
ax.XLim = [138 152]; % NE Japan
ax.YLim = [ 32  48]; % NE Japan
% ax.XLim = [130 142]; % SW Japan
% ax.YLim = [ 30  37]; % SW Japan
hold on

% Block line
for nb = 1:blk(1).nblock
  plot(blk(nb).lon,blk(nb).lat, 'm'); hold on
end

%% Pick edge line of locking segments
alon = [];
alat = [];
count = 0;
while 1
  [x,y,id1] = ginput(1);
  if id1 == 1
    count = count + 1;
    xc = x;
    yc = y;
    alon = [alon; [count, xc]];
    alat = [alat; [count, yc]];
    plot(xc,yc,'ob')
    text(xc,yc,[num2str(count)])
  else
    break
  end
end

% Add start point
alon = [alon; alon(1,:)];
alat = [alat; alat(1,:)];

% Test plot
plot(alon(:,2),alat(:,2),'b')

% Save
fid = fopen('patch_output.txt','wt');
fprintf(fid,'%8.3f %7.3f\n',[alon(:,2),alat(:,2)]');
fclose(fid);