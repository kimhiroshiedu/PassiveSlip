function[dSsn,dSdn,dUxyz]=makeG_d_04(triC,tri3)

%change the value%
pr=0.25;
ss= 0.0;
ds= 1.0;
ts= 0.0;

%Define strike and dip of fault%
%ca=[500,0,0];
%cb=[0,300,0];
%normVec                      = -cross(ca,cb);
%normVec                      = normVec./norm(normVec)
%if (normVec(3) < 0) % Enforce clockwise circulation
%   normVec                   = -normVec;
%   [x(2) x(3)]               = swap(x(2), x(3));
%   [y(2) y(3)]               = swap(y(2), y(3));
%   [z(2) z(3)]               = swap(z(2), z(3));
%end
%strikeVec                    = [-sin(atan2(normVec(2),normVec(1))),cos(atan2(normVec(2),normVec(1))),0]
%dipVec                       = cross(normVec, strikeVec)

%if normVec(1)==0 && normVec(2)==0
    strikeVec=[1 0 0]
    dipVec   =[0 1 0]
    normVec  =[0 0 -1]
%end

%RsitaS=atan2(strikeVec(2),strikeVec(1)); %from x-axis countour clockwise [rad]%
%RsitaSY=2.5*pi-RsitaS; %from y-axis clockwise[rad]%
%sitaSY=RsitaSY*180/pi
%sitaS=RsitaS;

%Rdip=asin(dipVec(3));
%dip=Rdip*180.0/pi;
%sitaD=Rdip+pi;
  
  n=length(triC);

  sx=triC(:,1);
  sy=triC(:,2);
  sz=triC(:,3)+100.0;

%  c1=cos(sitaS);
%  s1=sin(sitaS);
%  c2=cos(sitaD);
%  s2=sin(sitaD);

  Usumx=zeros(n,n);
  Usumy=zeros(n,n);
  Usumz=zeros(n,n);
  dUxyz=zeros(n,n,3);
  %Sxx=zeros(n,n);
  %Sxy=zeros(n,n);
  %Sxz=zeros(n,n);
  %Syy=zeros(n,n);
  %Syz=zeros(n,n);
  %Szz=zeros(n,n);
  
  dSsn=zeros(n,n);
  dSdn=zeros(n,n);
  
 for i=1:n%numbers of fault%
  x=[tri3(i,1,1),tri3(i,1,2),tri3(i,1,3)];
  y=[tri3(i,2,1),tri3(i,2,2),tri3(i,2,3)];
  z=[tri3(i,3,1),tri3(i,3,2),tri3(i,3,3)];
  
  [U] = CalcTriDisps(sx, sy, sz, x, y, z, pr, ss, ts, ds);
  Usumx(:,i)=U.x;
  Usumy(:,i)=U.y;
  Usumz(:,i)=-U.z;
  dUxyz(:,i,1)=Usumx(:,i);
  dUxyz(:,i,2)=Usumy(:,i);
  dUxyz(:,i,3)=Usumz(:,i);
  
  [S] = CalcTriStrains(sx, sy, sz, x, y, z, pr, ss, ts, ds);
  %Sxx(:,i)=S.xx;
  %Sxy(:,i)=S.xy;
  dSsn(:,i)=S.xz;
  %Syy(:,i)=S.yy;
  dSdn(:,i)=S.yz;
  %Szz(:,i)=S.zz;
  
  %s=strike, d=dip, n=normal
  %Sss(:,i)=c1^2*Sxx+2*c1*s1*Sxy+s1^2*Syy;
  %Ssd(:,i)=c2*(c1*s1*(Syy-Sxx)+(c1^2-s1^2)*Sxy)+s2*(c1*Sxz+s1*Syz);
  %dSsn(:,i)=-s2*(c1*s1*(Syy(:,i)-Sxx(:,i))+(c1^2-s1^2)*Sxy(:,i))+c2*(c1*Sxz(:,i)+s1*Syz(:,i));
  %Sdd(:,i)=c2^2*(s1^2*Sxx-2*s1*c1*Sxy+c1^2*Syy)+2*c2*s2*(c1*Syz-s1*Sxz)+s2^2*Szz;
  %dSdn(:,i)=c2*s2*(Szz(:,i)-s1^2*Sxx(:,i)+2*s1*c1*Sxy(:,i)-c1^2*Syy(:,i))+(c2^2-s2^2)*(c1*Syz(:,i)-s1*Sxz(:,i));
  %Snn(:,i)=s2^2*(s1^2*Sxx-2*s1*c1*Sxy+c1^2*Syy)-2*c2*s2*(c1*Syz-s1*Sxz)+c2^2*Szz;
 end
end
