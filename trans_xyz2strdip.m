%% Transform tensor from xyz to strike-dip
function [Gs]=trans_xyz2strdip(Gs,sitaS,sitaD)
% This function transforms a strain tensor from xyz to fault strike-dip.
% 
% Output
% Gs.stst : strain of strike direction on the fault due to strike slip.
% Gs.stdp : strain of dip direction on the fault due to strike slip.
% Gs.stts : strain of tensile direction on the fault due to strike slip.
% Gs.dpst : strain of strike direction on the fault due to dip slip.
% Gs.dpdp : strain of dip direction on the fault due to dip slip.
% Gs.dpts : strain of tensile direction on the fault due to dip slip.
% Gs.tsst : strain of strike direction on the fault due to tensile slip.
% Gs.tsdp : strain of dip direction on the fault due to tensile slip.
% Gs.tsts : strain of tensile direction on the fault due to tensile slip.
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

[sst,sdp,sts]=calctrans(Gs.st,c1,s1,c2,s2); % response to strike slip
Gs.stst=sst;
Gs.stdp=sdp;
Gs.stts=sts;
[sst,sdp,sts]=calctrans(Gs.ts,c1,s1,c2,s2); % response to tensile slip
Gs.tsst=sst;
Gs.tsdp=sdp;
Gs.tsts=sts;
[sst,sdp,sts]=calctrans(Gs.dp,c1,s1,c2,s2); % response to dip slip
Gs.dpst=sst;
Gs.dpdp=sdp;
Gs.dpts=sts;

end