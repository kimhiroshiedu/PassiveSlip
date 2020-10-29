%% make Green's function of displacement
function [Gu]=makeGreenDisp(xyz,trixyz3)
% This function make Green's function of displacement (Gu).
% Input
%  xyz     : site location.
%  trixyz3 : trimesh coordianate.
% Output
%  Gu.st   : surface displacement due to strike slip on a fault.
%  Gu.dp   : surface displacement due to dip slip on a fault.
%  Gu.ts   : surface displacement due to tensile slip on a fault.
%  *Each matrix (Gu.*) contain 3 x NOBS by NFLT elements.

pr  =0.25; % Poisson's ratio
nobs=size(xyz,1);
nflt=size(trixyz3,1);
sx  =xyz(:,1);
sy  =xyz(:,2);
sz  =zeros(nobs,1);

for n=1:nflt
    trix=[trixyz3(n,1),trixyz3(n,4),trixyz3(n,7)];
    triy=[trixyz3(n,2),trixyz3(n,5),trixyz3(n,8)];
    triz=[trixyz3(n,3),trixyz3(n,6),trixyz3(n,9)];
%     [U] = CalcTriDisps(sx, sy, sz, x, y, z, pr, ss, ts, ds);
    [U] = CalcTriDisps(sx, sy, sz, trix, triy, triz, pr, 1, 0, 0); % Strike
    Gu.st(1:3:3*nobs,n)= U.x;
    Gu.st(2:3:3*nobs,n)= U.y;
    Gu.st(3:3:3*nobs,n)=-U.z;
    [U] = CalcTriDisps(sx, sy, sz, trix, triy, triz, pr, 0, 1, 0); % Tensile
    Gu.ts(1:3:3*nobs,n)= U.x;
    Gu.ts(2:3:3*nobs,n)= U.y;
    Gu.ts(3:3:3*nobs,n)=-U.z;
    [U] = CalcTriDisps(sx, sy, sz, trix, triy, triz, pr, 0, 0, 1); % Dip
    Gu.dp(1:3:3*nobs,n)= U.x;
    Gu.dp(2:3:3*nobs,n)= U.y;
    Gu.dp(3:3:3*nobs,n)=-U.z;
end

end
