function pick_locksegments_from_result(blk,tcha,d,G,T)
% before using this script, load blk, tcha, grn,
G(1).tb_kin = full(G(1).tb_kin);
G(1).tb_mec = full(G(1).tb_mec);

% Calculate back-slip on locked and creeping patches.
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

% Make figure showing mesh, coast and block lines, and slip rate
% Additionally, calculate the center and sdr on the mesh
figure(100); clf(100)
clon = [];
clat = [];
cdep = [];
csdr = [];
crel = [];
mm3 = 1;
for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    if blk(1).bound(nb1,nb2).flag2 == 1
      nf = size(blk(1).bound(nb1,nb2).blon,1);
      % Backslip rate (locked)
      patch(blk(1).bound(nb1,nb2).blon',blk(1).bound(nb1,nb2).blat',...
            single(sqrt(bslipl50(mm3:mm3+nf-1,1).^2+bslipl50(mm3+nf:mm3+2*nf-1,1).^2)~=0));
      % Backslip rate (total)
      %       patch(blk(1).bound(nb1,nb2).blon',blk(1).bound(nb1,nb2).blat',...
      %             sqrt(bslip50(mm3:mm3+nf-1,1).^2+bslip50(mm3+nf:mm3+2*nf-1,1).^2));
      hold on
      clon = [clon; mean(blk(1).bound(nb1,nb2).blon,2)];
      clat = [clat; mean(blk(1).bound(nb1,nb2).blat,2)];
      cdep = [cdep; mean(blk(1).bound(nb1,nb2).bdep,2)];
      csdr = [csdr; sqrt(bslip(mm3:mm3+nf-1,1).^2+bslip(mm3+nf:mm3+2*nf-1,1).^2)];
      crel = [crel; sqrt(rel_mec(mm3:mm3+nf-1,1).^2+rel_mec(mm3+nf:mm3+2*nf-1,1).^2)];
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
% ax.XLim = [138 152]; % NE Japan
% ax.YLim = [ 32  48]; % NE Japan
ax.XLim = [130 142]; % SW Japan
ax.YLim = [ 30  37]; % SW Japan
hold on

% Block line
for nb = 1:blk(1).nblock
  plot(blk(nb).lon,blk(nb).lat, 'm'); hold on
end

end

%% Pick edge line of locking segments
function pick_points
alon = [];
alat = [];
count = 0;
while 1
  [x_local,y_local,id1] = ginput(1);
  if id1 == 1
    count = count + 1;
    xc = x_local;
    yc = y_local;
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

end

%%
function calc_sdrline(clon,clat,cdep,csdr,crel)
% Prepare calc line of SDR and Vplate
alon2 = linspace2(alon(:,2),100);
alat2 = linspace2(alat(:,2),100);
[x_local,y_local] = PLTXY(alon2,alat2,alon2(end),alat2(end));
leng = sqrt((x_local-x_local(end)).^2 + (y_local-y_local(end)).^2);
F1 = scatteredInterpolant(clon,clat,double(csdr));
F2 = scatteredInterpolant(clon,clat,double(crel));
F3 = scatteredInterpolant(clon,clat,cdep);
asdr2 = F1(alon2,alat2);
arel2 = F2(alon2,alat2);
adep2 = F3(alon2,alat2);
fid = fopen('slip_rate_along_line.txt','wt');
fprintf(fid,'# Lon  Lat  Dep  Leng  SDR  Vpl\n');
fprintf(fid,'%8.3f %7.3f %6.2f %9.3f %6.2f %6.2f\n', [alon2,alat2,adep2,leng,asdr2,arel2]');
fclose(fid);

figure(200); clf(200);
p2 = plot(leng,arel2,"LineStyle","--","Color",'r');
hold on
p1 = plot(leng,asdr2,"LineStyle","-","Color",'k',"LineWidth",1.5);
ax = gca;
ax.FontName = 'Helvetica';
ax.FontSize = 14;
ax.XLabel.String = 'Distance (km)';
ax.YLabel.String = 'Rate (mm/yr)';
legend([p1 p2],{'Slip deficit','Plate convergence'})
ax.PlotBoxAspectRatio = [1.7,1,1];
print('slip_rate_along_line','-dpdf','-painters')

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

%% PLTXY
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
