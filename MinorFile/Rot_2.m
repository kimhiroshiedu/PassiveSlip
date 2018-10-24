function [R2Sxz,R2Syz]=Rot_2(sitaS,sitaD,Sxx,Sxy,Sxz,Syy,Syz,Szz)
c1=cos(sitaS);
s1=sin(sitaS);
c2=cos(sitaD);
s2=sin(sitaD);

%R2Sxx=c1^2*Sxx+2*c1*s1*Sxy+s1^2*Syy;
%R2Sxy=c2*(c1*s1*(Syy-Sxx)+(c1^2-s1^2)*Sxy)+s2*(c1*Sxz+s1*Syz);
R2Sxz=-s2*(c1*s1*(Syy-Sxx)+(c1^2-s1^2)*Sxy)+c2*(c1*Sxz+s1*Syz);
%R2Syy=c2^2*(s1^2*Sxx-2*s1*c1*Sxy+c1^2*Syy)+2*c2*s2*(c1*Syz-s1*Sxz)+s2^2*Szz;
R2Syz=c2*s2*(Szz-s1^2*Sxx+2*s1*c1*Sxy-c1^2*Syy)+(c2^2-s2^2)*(c1*Syz-s1*Sxz);
%R2Szz=s2^2*(s1^2*Sxx-2*s1*c1*Sxy+c1^2*Syy)-2*c2*s2*(c1*Syz-s1*Sxz)+c2^2*Szz;
end
