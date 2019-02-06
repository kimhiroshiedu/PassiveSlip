function [v1, v2, v3] = adv(y1, y2, y3, a, beta, nu, B1, B2, B3)
% These are the displacements in a uniform elastic half space due to slip
% on an angular dislocation (Comninou and Dunders, 1975).  Some of the
% equations for the B2 and B3 cases have been corrected following Thomas
% 1993.  The equations are coded in way such that they roughly correspond
% to each line in original text.  Exceptions have been made where it made
% more sense because of grouping symbols.
sinbeta           = sin(beta);
cosbeta           = cos(beta);
cotbeta           = cot(beta);
z1                = y1.*cosbeta - y3.*sinbeta;
z3                = y1.*sinbeta + y3.*cosbeta;
R2                = y1.*y1 + y2.*y2 + y3.*y3;
R                 = sqrt(R2);
y3bar             = y3 + 2.*a;
z1bar             = y1.*cosbeta + y3bar.*sinbeta;
z3bar             = -y1.*sinbeta + y3bar.*cosbeta;
R2bar             = y1.*y1 + y2.*y2 + y3bar.*y3bar;
Rbar              = sqrt(R2bar);
F                 = -atan2(y2, y1) + atan2(y2, z1) + atan2(y2.*R.*sinbeta, y1.*z1+(y2.*y2).*cosbeta);
Fbar              = -atan2(y2, y1) + atan2(y2, z1bar) + atan2(y2.*Rbar.*sinbeta, y1.*z1bar+(y2.*y2).*cosbeta);

% Case I: Burgers vector (B1,0,0)
v1InfB1           = 2.*(1-nu).*(F+Fbar) - y1.*y2.*(1./(R.*(R-y3)) + 1./(Rbar.*(Rbar+y3bar))) - ...
    y2.*cosbeta.*((R.*sinbeta-y1)./(R.*(R-z3)) + (Rbar.*sinbeta-y1)./(Rbar.*(Rbar+z3bar)));
v2InfB1           = (1-2.*nu).*(log(R-y3)+log(Rbar+y3bar) - cosbeta.*(log(R-z3)+log(Rbar+z3bar))) - ...
    y2.*y2.*(1./(R.*(R-y3))+1./(Rbar.*(Rbar+y3bar)) - cosbeta.*(1./(R.*(R-z3))+1./(Rbar.*(Rbar+z3bar))));
v3InfB1           = y2 .* (1./R - 1./Rbar - cosbeta.*((R.*cosbeta-y3)./(R.*(R-z3)) - (Rbar.*cosbeta+y3bar)./(Rbar.*(Rbar+z3bar))));
v1InfB1           = v1InfB1 ./ (8.*pi.*(1-nu));
v2InfB1           = v2InfB1 ./ (8.*pi.*(1-nu));
v3InfB1           = v3InfB1 ./ (8.*pi.*(1-nu));

v1CB1             = -2.*(1-nu).*(1-2.*nu).*Fbar.*(cotbeta.*cotbeta) + (1-2.*nu).*y2./(Rbar+y3bar) .* ((1-2.*nu-a./Rbar).*cotbeta - y1./(Rbar+y3bar).*(nu+a./Rbar)) + ...
    (1-2.*nu).*y2.*cosbeta.*cotbeta./(Rbar+z3bar).*(cosbeta+a./Rbar) + a.*y2.*(y3bar-a).*cotbeta./(Rbar.*Rbar.*Rbar) + ...
    y2.*(y3bar-a)./(Rbar.*(Rbar+y3bar)).*(-(1-2.*nu).*cotbeta + y1./(Rbar+y3bar) .* (2.*nu+a./Rbar) + a.*y1./(Rbar.*Rbar)) + ...
    y2.*(y3bar-a)./(Rbar.*(Rbar+z3bar)).*(cosbeta./(Rbar+z3bar).*((Rbar.*cosbeta+y3bar) .* ((1-2.*nu).*cosbeta-a./Rbar).*cotbeta + 2.*(1-nu).*(Rbar.*sinbeta-y1).*cosbeta) - a.*y3bar.*cosbeta.*cotbeta./(Rbar.*Rbar));
v2CB1             = (1-2.*nu).*((2.*(1-nu).*(cotbeta.*cotbeta)-nu).*log(Rbar+y3bar) -(2.*(1-nu).*(cotbeta.*cotbeta)+1-2.*nu).*cosbeta.*log(Rbar+z3bar)) - ...
    (1-2.*nu)./(Rbar+y3bar).*(y1.*cotbeta.*(1-2.*nu-a./Rbar) + nu.*y3bar - a + (y2.*y2)./(Rbar+y3bar).*(nu+a./Rbar)) - ...
    (1-2.*nu).*z1bar.*cotbeta./(Rbar+z3bar).*(cosbeta+a./Rbar) - a.*y1.*(y3bar-a).*cotbeta./(Rbar.*Rbar.*Rbar) + ...
    (y3bar-a)./(Rbar+y3bar).*(-2.*nu + 1./Rbar.*((1-2.*nu).*y1.*cotbeta-a) + (y2.*y2)./(Rbar.*(Rbar+y3bar)).*(2.*nu+a./Rbar)+a.*(y2.*y2)./(Rbar.*Rbar.*Rbar)) + ...
    (y3bar-a)./(Rbar+z3bar).*((cosbeta.*cosbeta) - 1./Rbar.*((1-2.*nu).*z1bar.*cotbeta+a.*cosbeta) + a.*y3bar.*z1bar.*cotbeta./(Rbar.*Rbar.*Rbar) - 1./(Rbar.*(Rbar+z3bar)) .* ((y2.*y2).*(cosbeta.*cosbeta) - a.*z1bar.*cotbeta./Rbar.*(Rbar.*cosbeta+y3bar)));

v3CB1             = 2.*(1-nu).*(((1-2.*nu).*Fbar.*cotbeta) + (y2./(Rbar+y3bar).*(2.*nu+a./Rbar)) - (y2.*cosbeta./(Rbar+z3bar).*(cosbeta+a./Rbar))) + ...
    y2.*(y3bar-a)./Rbar.*(2.*nu./(Rbar+y3bar)+a./(Rbar.*Rbar)) + ...
    y2.*(y3bar-a).*cosbeta./(Rbar.*(Rbar+z3bar)).*(1-2.*nu-(Rbar.*cosbeta+y3bar)./(Rbar+z3bar).*(cosbeta + a./Rbar) - a.*y3bar./(Rbar.*Rbar));

v1CB1             = v1CB1 ./ (4.*pi.*(1-nu));
v2CB1             = v2CB1 ./ (4.*pi.*(1-nu));
v3CB1             = v3CB1 ./ (4.*pi.*(1-nu));

v1B1              = v1InfB1 + v1CB1;
v2B1              = v2InfB1 + v2CB1;
v3B1              = v3InfB1 + v3CB1;


% Case II: Burgers vector (0,B2,0)
v1InfB2           = -(1-2.*nu).*(log(R-y3) + log(Rbar+y3bar)-cosbeta.*(log(R-z3)+log(Rbar+z3bar))) + ...
    y1.*y1.*(1./(R.*(R-y3))+1./(Rbar.*(Rbar+y3bar))) + z1.*(R.*sinbeta-y1)./(R.*(R-z3)) + z1bar.*(Rbar.*sinbeta-y1)./(Rbar.*(Rbar+z3bar));
v2InfB2           = 2.*(1-nu).*(F+Fbar) + y1.*y2.*(1./(R.*(R-y3))+1./(Rbar.*(Rbar+y3bar))) - y2.*(z1./(R.*(R-z3))+z1bar./(Rbar.*(Rbar+z3bar)));
v3InfB2           = -(1-2.*nu).*sinbeta.*(log(R-z3)-log(Rbar+z3bar)) - y1.*(1./R-1./Rbar) + z1.*(R.*cosbeta-y3)./(R.*(R-z3)) - z1bar.*(Rbar.*cosbeta+y3bar)./(Rbar.*(Rbar+z3bar));
v1InfB2           = v1InfB2 ./ (8.*pi.*(1-nu));
v2InfB2           = v2InfB2 ./ (8.*pi.*(1-nu));
v3InfB2           = v3InfB2 ./ (8.*pi.*(1-nu));

v1CB2             = (1-2.*nu).*((2.*(1-nu).*(cotbeta.*cotbeta)+nu).*log(Rbar+y3bar) - (2.*(1-nu).*(cotbeta.*cotbeta)+1).*cosbeta.*log(Rbar+z3bar)) + ...
    (1-2.*nu)./(Rbar+y3bar).* (-(1-2.*nu).*y1.*cotbeta+nu.*y3bar-a+a.*y1.*cotbeta./Rbar + (y1.*y1)./(Rbar+y3bar).*(nu+a./Rbar)) - ...
    (1-2.*nu).*cotbeta./(Rbar+z3bar).*(z1bar.*cosbeta - a.*(Rbar.*sinbeta-y1)./(Rbar.*cosbeta)) - a.*y1.*(y3bar-a).*cotbeta./(Rbar.*Rbar.*Rbar) + ...
    (y3bar-a)./(Rbar+y3bar).*(2.*nu + 1./Rbar.*((1-2.*nu).*y1.*cotbeta+a) - (y1.*y1)./(Rbar.*(Rbar+y3bar)).*(2.*nu+a./Rbar) - a.*(y1.*y1)./(Rbar.*Rbar.*Rbar)) + ...
    (y3bar-a).*cotbeta./(Rbar+z3bar).*(-cosbeta.*sinbeta+a.*y1.*y3bar./(Rbar.*Rbar.*Rbar.*cosbeta) + (Rbar.*sinbeta-y1)./Rbar.*(2.*(1-nu).*cosbeta - (Rbar.*cosbeta+y3bar)./(Rbar+z3bar).*(1+a./(Rbar.*cosbeta))));
v2CB2             = 2.*(1-nu).*(1-2.*nu).*Fbar.*cotbeta.*cotbeta + (1-2.*nu).*y2./(Rbar+y3bar).*(-(1-2.*nu-a./Rbar).*cotbeta + y1./(Rbar+y3bar).*(nu+a./Rbar)) - ...
    (1-2.*nu).*y2.*cotbeta./(Rbar+z3bar).*(1+a./(Rbar.*cosbeta)) - a.*y2.*(y3bar-a).*cotbeta./(Rbar.*Rbar.*Rbar) + ...
    y2.*(y3bar-a)./(Rbar.*(Rbar+y3bar)).*((1-2.*nu).*cotbeta - 2.*nu.*y1./(Rbar+y3bar) - a.*y1./Rbar.*(1./Rbar+1./(Rbar+y3bar))) + ...
    y2.*(y3bar-a).*cotbeta./(Rbar.*(Rbar+z3bar)).*(-2.*(1-nu).*cosbeta + (Rbar.*cosbeta+y3bar)./(Rbar+z3bar).*(1+a./(Rbar.*cosbeta)) + a.*y3bar./((Rbar.*Rbar).*cosbeta));
v3CB2             = -2.*(1-nu).*(1-2.*nu).*cotbeta .* (log(Rbar+y3bar)-cosbeta.*log(Rbar+z3bar)) - ...
    2.*(1-nu).*y1./(Rbar+y3bar).*(2.*nu+a./Rbar) + 2.*(1-nu).*z1bar./(Rbar+z3bar).*(cosbeta+a./Rbar) + ...
    (y3bar-a)./Rbar.*((1-2.*nu).*cotbeta-2.*nu.*y1./(Rbar+y3bar)-a.*y1./(Rbar.*Rbar)) - ...
    (y3bar-a)./(Rbar+z3bar).*(cosbeta.*sinbeta + (Rbar.*cosbeta+y3bar).*cotbeta./Rbar.*(2.*(1-nu).*cosbeta - (Rbar.*cosbeta+y3bar)./(Rbar+z3bar)) + a./Rbar.*(sinbeta - y3bar.*z1bar./(Rbar.*Rbar) - z1bar.*(Rbar.*cosbeta+y3bar)./(Rbar.*(Rbar+z3bar))));
v1CB2             = v1CB2 ./ (4.*pi.*(1-nu));
v2CB2             = v2CB2 ./ (4.*pi.*(1-nu));
v3CB2             = v3CB2 ./ (4.*pi.*(1-nu));

v1B2              = v1InfB2 + v1CB2;
v2B2              = v2InfB2 + v2CB2;
v3B2              = v3InfB2 + v3CB2;


% Case III: Burgers vector (0,0,B3)
v1InfB3           = y2.*sinbeta.*((R.*sinbeta-y1)./(R.*(R-z3))+(Rbar.*sinbeta-y1)./(Rbar.*(Rbar+z3bar)));
v2InfB3           = (1-2.*nu).*sinbeta.*(log(R-z3)+log(Rbar+z3bar)) - (y2.*y2).*sinbeta.*(1./(R.*(R-z3))+1./(Rbar.*(Rbar+z3bar)));
v3InfB3           = 2.*(1-nu).*(F-Fbar) + y2.*sinbeta.*((R.*cosbeta-y3)./(R.*(R-z3))-(Rbar.*cosbeta+y3bar)./(Rbar.*(Rbar+z3bar)));
v1InfB3           = v1InfB3 ./ (8.*pi.*(1-nu));
v2InfB3           = v2InfB3 ./ (8.*pi.*(1-nu));
v3InfB3           = v3InfB3 ./ (8.*pi.*(1-nu));

v1CB3             = (1-2.*nu).*(y2./(Rbar+y3bar).*(1+a./Rbar) - y2.*cosbeta./(Rbar+z3bar).*(cosbeta+a./Rbar)) - ...
    y2.*(y3bar-a)./Rbar.*(a./(Rbar.*Rbar) + 1./(Rbar+y3bar)) + ...
    y2.*(y3bar-a).*cosbeta./(Rbar.*(Rbar+z3bar)).*((Rbar.*cosbeta+y3bar)./(Rbar+z3bar).*(cosbeta+a./Rbar) + a.*y3bar./(Rbar.*Rbar));
v2CB3             = (1-2.*nu).*(-sinbeta.*log(Rbar+z3bar) - y1./(Rbar+y3bar).*(1+a./Rbar) + z1bar./(Rbar+z3bar).*(cosbeta+a./Rbar)) + ...
    y1.*(y3bar-a)./Rbar.*(a./(Rbar.*Rbar) + 1./(Rbar+y3bar)) - ...
    (y3bar-a)./(Rbar+z3bar).*(sinbeta.*(cosbeta-a./Rbar) + z1bar./Rbar.*(1+a.*y3bar./(Rbar.*Rbar)) - ...
    1./(Rbar.*(Rbar+z3bar)).*((y2.*y2).*cosbeta.*sinbeta - a.*z1bar./Rbar.*(Rbar.*cosbeta+y3bar)));
v3CB3             = 2.*(1-nu).*Fbar + 2.*(1-nu).*(y2.*sinbeta./(Rbar+z3bar).*(cosbeta + a./Rbar)) + ...
    y2.*(y3bar-a).*sinbeta./(Rbar.*(Rbar+z3bar)).*(1 + (Rbar.*cosbeta+y3bar)./(Rbar+z3bar).*(cosbeta+a./Rbar) + a.*y3bar./(Rbar.*Rbar));
v1CB3             = v1CB3 ./ (4.*pi.*(1-nu));
v2CB3             = v2CB3 ./ (4.*pi.*(1-nu));
v3CB3             = v3CB3 ./ (4.*pi.*(1-nu));

v1B3              = v1InfB3 + v1CB3;
v2B3              = v2InfB3 + v2CB3;
v3B3              = v3InfB3 + v3CB3;


% Sum the for each slip component
v1                = B1.*v1B1 + B2.*v1B2 + B3.*v1B3;
v2                = B1.*v2B1 + B2.*v2B2 + B3.*v2B3;
v3                = B1.*v3B1 + B2.*v3B2 + B3.*v3B3;
end