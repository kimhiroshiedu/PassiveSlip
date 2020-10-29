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
[Gs]=trans_xyz2strdip(Gs,sitaS,sitaD);
end