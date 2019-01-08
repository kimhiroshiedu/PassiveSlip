%% Transform tensor from xyz to strike-dip
function[gsf]=trans_xyz2strdip(gsx,sitaS,sitaD)
% This function transforms a strain tensor from xyz to fault strike-dip.
% 
% E_stdp = Tx * Tz * E_xyz * Tz' * Tx'
% <->
% | exx' exy' exz' |   | 1   0  0 |   |  c1 s1 0 |   | exx exy exz |   |  c1 s1 0 |T  | 1   0  0 |T
% | eyx' eyy' eyz' | = | 0  c2 s2 | * | -s1 c1 0 | * | eyx eyy eyz | * | -s1 c1 0 | * | 0  c2 s2 |
% | ezx' ezy' ezz' |   | 0 -s2 c2 |   |   0  0 1 |   | ezx ezy ezz |   |   0  0 1 |   | 0 -s2 c2 |
% where 
% c1 : cos(strike)
% s1 : sin(strike)
% c2 : cos(dip)
% s2 : sin(dip)
% Note strike is the counter clock wise angle from x-axis, dip is angle
% from xy-plane. "T" indicates the transpose of matrix.
% 
% Output
% gsf.ss : strain of strike direction on the fault due to strike slip.
% gsf.sd : strain of dip direction on the fault due to strike slip.
% gsf.st : strain of tensile direction on the fault due to strike slip.
% gsf.ds : strain of strike direction on the fault due to dip slip.
% gsf.dd : strain of dip direction on the fault due to dip slip.
% gsf.dt : strain of tensile direction on the fault due to dip slip.
% gsf.ts : strain of strike direction on the fault due to tensile slip.
% gsf.td : strain of dip direction on the fault due to tensile slip.
% gsf.tt : strain of tensile direction on the fault due to tensile slip.
% 
% Note that strain of strike direction on the fault corresponds to ezx',
% strain of dip direction on the fault corresponds to ezy', and
% strain of tensile direction on the fault corresponds to ezz'.

% c1=cos(sitaS);
% s1=sin(sitaS);
% c2=cos(sitaD);
% s2=sin(sitaD);
c1=repmat(cos(sitaS)',1,size(sitaS,2));
s1=repmat(sin(sitaS)',1,size(sitaS,2));
c2=repmat(cos(sitaD)',1,size(sitaD,2));
s2=repmat(sin(sitaD)',1,size(sitaD,2));

[sst,sdp,sts]=calctrans(gsx.st,c1,s1,c2,s2); % response to strike slip
gsf.ss=sst;
gsf.sd=sdp;
gsf.st=sts;
[sst,sdp,sts]=calctrans(gsx.ts,c1,s1,c2,s2); % response to tensile slip
gsf.ts=sst;
gsf.td=sdp;
gsf.tt=sts;
[sst,sdp,sts]=calctrans(gsx.dp,c1,s1,c2,s2); % response to dip slip
gsf.ds=sst;
gsf.dd=sdp;
gsf.dt=sts;

end
