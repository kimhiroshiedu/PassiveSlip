%% make Green's function of strain
function[gs]=makeGreenStrain(trixyzC,trixyz3,sitaS,sitaD,normVec)
% This function make Green's function of strain (gs).
% Input
%  trixyzC : center of trimesh.
%  trixyz3 : trimesh coordianate.
%  sitaS   : angle of fault strike from X-axis.
%  sitaD   : angle of fault dip from XY-plane.
%  normVec : normalized normal vectors of faults.
% Output
%  gs.st   : strain due to strike slip on a fault.
%  gs.dp   : strain due to dip slip on a fault.
%  gs.ts   : strain due to tensile slip on a fault.
%  *Each matrix (gs.*) contain 6 x NFLT by NFLT elements.
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
    gs.st(1:6:6*nflt,n)=S.xx;
    gs.st(2:6:6*nflt,n)=S.xy;
    gs.st(3:6:6*nflt,n)=S.xz;
    gs.st(4:6:6*nflt,n)=S.yy;
    gs.st(5:6:6*nflt,n)=S.yz;
    gs.st(6:6:6*nflt,n)=S.zz;
    [S] = CalcTriStrains(sx, sy, sz, trix, triy, triz, pr, 0, 1, 0); % Tensile
    gs.ts(1:6:6*nflt,n)=S.xx;
    gs.ts(2:6:6*nflt,n)=S.xy;
    gs.ts(3:6:6*nflt,n)=S.xz;
    gs.ts(4:6:6*nflt,n)=S.yy;
    gs.ts(5:6:6*nflt,n)=S.yz;
    gs.ts(6:6:6*nflt,n)=S.zz;
    [S] = CalcTriStrains(sx, sy, sz, trix, triy, triz, pr, 0, 0, 1); % Dip
    gs.dp(1:6:6*nflt,n)=S.xx;
    gs.dp(2:6:6*nflt,n)=S.xy;
    gs.dp(3:6:6*nflt,n)=S.xz;
    gs.dp(4:6:6*nflt,n)=S.yy;
    gs.dp(5:6:6*nflt,n)=S.yz;
    gs.dp(6:6:6*nflt,n)=S.zz;
end

end
