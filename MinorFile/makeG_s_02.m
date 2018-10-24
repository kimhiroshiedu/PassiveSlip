function[sR2Sxz,sR2Syz]=makeG_s_02(triC,tri3)

%change the value%
pr=0.25;
ss= 1;
ds= 0;
ts= 0;

%Define strike and dip of fault%
ca=[500,0,0];
cb=[0,300,300*-0.3];
normVec                      = cross(ca,cb);
normVec                      = normVec./norm(normVec)
%if (normVec(3) < 0) % Enforce clockwise circulation
%   normVec                   = -normVec;
%   [x(2) x(3)]               = swap(x(2), x(3));
%   [y(2) y(3)]               = swap(y(2), y(3));
%   [z(2) z(3)]               = swap(z(2), z(3));
%end
strikeVec                    = [-sin(atan2(normVec(2),normVec(1))),cos(atan2(normVec(2),normVec(1))),0]
dipVec                       = -cross(normVec, strikeVec)

%if normVec(1)==0 and normVec(2)==0
%    strikeVec=[1 0 0]
%    dipVec   =[0 1 0]
%    normVec  =[0 0 1]
%end

RsitaS=atan2(strikeVec(2),strikeVec(1)); %from x-axis countour clockwise [rad]%
RsitaSY=2.5*pi-RsitaS; %from y-axis clockwise[rad]%
sitaSY=RsitaSY*180/pi
sitaS=RsitaS+pi;

Rdip=asin(-dipVec(3));
dip=Rdip*180/pi
sitaD=2.0*pi-Rdip;
  
  n=length(triC);

  sx=triC(:,1);
  sy=triC(:,2);
  sz=triC(:,3)-500;

  c1=cos(sitaS);
  s1=sin(sitaS);
  c2=cos(sitaD);
  s2=sin(sitaD);

  
 for i=1:n%numbers of fault%
  x=[tri3(i,1,1),tri3(i,1,2),tri3(i,1,3)];
  y=[tri3(i,2,1),tri3(i,2,2),tri3(i,2,3)];
  z=[tri3(i,3,1),tri3(i,3,2),tri3(i,3,3)];
  
  %[U] = CalcTriDisps(sx, sy, sz, x, y, z, pr, ss, ts, ds);
  %Usumx=U.x;
  %Usumy=U.y;
  %Usumz=-U.z;
  [S] = CalcTriStrains(sx, sy, sz, x, y, z, pr, ss, ts, ds);
  Sxx(:,i)=S.xx;
  Sxy(:,i)=S.xy;
  Sxz(:,i)=S.xz;
  Syy(:,i)=S.yy;
  Syz(:,i)=S.yz;
  Szz(:,i)=S.zz;

  %R2Sxx(:,i)=c1^2*Sxx+2*c1*s1*Sxy+s1^2*Syy;
  %R2Sxy(:,i)=c2*(c1*s1*(Syy-Sxx)+(c1^2-s1^2)*Sxy)+s2*(c1*Sxz+s1*Syz);
  sR2Sxz(:,i)=-s2*(c1*s1*(Syy(:,i)-Sxx(:,i))+(c1^2-s1^2)*Sxy(:,i))+c2*(c1*Sxz(:,i)+s1*Syz(:,i));
  %R2Syy(:,i)=c2^2*(s1^2*Sxx-2*s1*c1*Sxy+c1^2*Syy)+2*c2*s2*(c1*Syz-s1*Sxz)+s2^2*Szz;
  sR2Syz(:,i)=c2*s2*(Szz(:,i)-s1^2*Sxx(:,i)+2*s1*c1*Sxy(:,i)-c1^2*Syy(:,i))+(c2^2-s2^2)*(c1*Syz(:,i)-s1*Sxz(:,i));
  %R2Szz(:,i)=s2^2*(s1^2*Sxx-2*s1*c1*Sxy+c1^2*Syy)-2*c2*s2*(c1*Syz-s1*Sxz)+c2^2*Szz;
 end 
end
