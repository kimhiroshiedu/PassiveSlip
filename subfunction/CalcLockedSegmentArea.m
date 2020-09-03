function CalcLockedSegmentArea(savefolder,blk,obs,tcha,patchfolder)
blk = ReadLockedPatch(blk,patchfolder);
SaveAsperitySegmentArea(savefolder,blk,obs,tcha);
end

%%
function SaveAsperitySegmentArea(folder,blk,obs,tcha)
savedir = fullfile(folder,'backslip');
if exist(savedir) ~=7; mkdir(savedir); end
alat0 = mean(obs(1).alat,2);
alon0 = mean(obs(1).alon,2);
mm1m = 1;
mm3m = 1;
mm1k = 1;
mm3k = 1;
mm1  = 1;
mm3  = 1;
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
      blk(1).bound(nb1,nb2).triarea = area;
      if blk(1).bound(nb1,nb2).flag2 == 1
        if blk(1).bound(nb1,nb2).segmentid == 1
          for nseg = 1:size(blk(1).bound(nb1,nb2).segment,2)
            [edgex,edgey] = PLTXY(blk(1).bound(nb1,nb2).segment(nseg).lat,blk(1).bound(nb1,nb2).segment(nseg).lon,alat0,alon0);
            segid = inpolygon(tricx,tricy,edgex,edgey);
            asp(1).bound(nb1,nb2).segment(nseg).smparea = segid'.*blk(1).bound(nb1,nb2).triarea' * tcha.smpaid(mm1m:mm1m+nf-1,:);
          end
        else
          asp(1).bound(nb1,nb2).smparea = blk(1).bound(nb1,nb2).triarea' * tcha.smpaid(mm1m:mm1m+nf-1,:);
        end
        mm1m = mm1m +   nf;
        mm3m = mm3m + 3*nf;  
      else
        mm1k = mm1k +   nf;
        mm3k = mm3k + 3*nf;          
      end
      mm1 = mm1 +   nf;
      mm3 = mm3 + 3*nf;
    end
  end
end
save(fullfile(savedir,'asp'),'asp');
end

%% Read locked patches
function [blk] = ReadLockedPatch(blk,patchfolder)
%    Test version coded by H. Kimura 2019/1/29
% Revised version coded by H. Kimura 2019/2/5

for nb1 = 1:blk(1).nblock
  for nb2 = nb1+1:blk(1).nblock
    if blk(1).bound(nb1,nb2).flag2 == 1
      blk(1).bound(nb1,nb2).segmentid = 0;
      patchfile = fullfile(patchfolder,['segments_',num2str(nb1),'_',num2str(nb2),'.txt']);
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
            blk(1).bound(nb1,nb2).segment(np+1).lon(n) = str2double(cellstr(tmp(1)));
            blk(1).bound(nb1,nb2).segment(np+1).lat(n) = str2double(cellstr(tmp(2)));
          else
            np = np+1;
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

%% Calculate triangle area
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