function SasaTriDisloca2
pi=3.14159265358979;
%change the value%
pr=0.25;
ss= 0;
ds= 0.001;
ts= 0;
%x=[xa,xb,xc],y=[ya,yb,yc],z=[za,zb,zc] vertices of triangle%
x=[-25,25,0];
y=[-25,-25,25];
z=[-35,-35,-50];
%make observe point on plane parallel to the triangle%
ca=[x(1,1)-x(1,3),y(1,1)-y(1,3),z(1,1)-z(1,3)];
cb=[x(1,2)-x(1,3),y(1,2)-y(1,3),z(1,2)-z(1,3)];
oc=[x(1,3),y(1,3),z(1,3)];

%Define strike and dip of fault%
normVec                      = cross(ca,cb);
normVec                      = normVec./norm(normVec)
if (normVec(3) < 0) % Enforce clockwise circulation
   normVec                   = -normVec;
   [x(2) x(3)]               = swap(x(2), x(3));
   [y(2) y(3)]               = swap(y(2), y(3));
   [z(2) z(3)]               = swap(z(2), z(3));
end
strikeVec                    = [-sin(atan2(normVec(2),normVec(1))),cos(atan2(normVec(2),normVec(1))),0]
dipVec                       = cross(normVec, strikeVec)

%if normVec(1)==0 and normVec(2)==0
%    strikeVec=[1 0 0]
%    dipVec   =[0 1 0]
%    normVec  =[0 0 1]
%end

RsitaS=atan2(strikeVec(2),strikeVec(1)); %from x-axis countour clockwise [rad]%
RsitaSY=2.5*pi-RsitaS; %from y-axis clockwise[rad]%
sitaSY=RsitaSY*180/pi

Rdip=asin(dipVec(3));
dip=Rdip*180/pi


%make Roatation Matrix%
RsitaSX90=RsitaS+pi; %countour clockwise from x-axis +pi [rad]%

RotS=[cos(RsitaSX90),-sin(RsitaSX90),0;sin(RsitaSX90),cos(RsitaSX90),0;0,0,1];

Rdip360=2.0*pi-Rdip; %countour clockwise from earth surface%

RotD=[1,0,0;0,cos(Rdip360),-sin(Rdip360);0,sin(Rdip360),cos(Rdip360)];

%make observed point as parallelogram%
i=0;

pp=-1/100;
 for p=1:100
  pp=pp+1/100;
  qq=-1/100;
   for q=1:100
    i=i+1;
    qq=qq+1/100;

    cp=[pp*ca(1,1)+qq*cb(1,1),pp*ca(1,2)+qq*cb(1,2),pp*ca(1,3)+qq*cb(1,3)];

    op=[oc(1,1)+cp(1,1),oc(1,2)+cp(1,2),oc(1,3)+cp(1,3)];

    sxyzplus(i,:)=[op(1,1),op(1,2),op(1,3)+0.01];
    %sxyzminus(i,:)=[op(1,1),op(1,2),op(1,3)-0.01];

   end
 end
%writing data%
%Fid1=fopen('/home_tmp/sasajima/DATA/001xydisp.dat','w');
%Fid2=fopen('/home/sasajima/DATA/005xystrain.dat','w');

%i=0;
%Calculation%
% for n=1:10000
%  i=i+1;
  sx=sxyzplus(:,1);
  sy=sxyzplus(:,2);
  sz=sxyzplus(:,3);
  %[U] = CalcTriDisps(sx, sy, sz, x, y, z, pr, ss, ts, ds);
  %Usumx=U.x;
  %Usumy=U.y;
  %Usumz=-U.z;
  [S] = CalcTriStrains(sx, sy, sz, x, y, z, pr, ss, ts, ds);
  Sxx=S.xx;
  Sxy=S.xy;
  Sxz=S.xz;
  Syy=S.yy;
  Syz=S.yz;
  Szz=S.zz;
for i=1:length(sx)
  Strain=[Sxx(i),Sxy(i),Sxz(i);Sxy(i),Syy(i),Syz(i);Sxz(i),Syz(i),Szz(i)];
  %StrainF=RotD'*RotS'*Strain*RotS*RotD;
  Strain1=RotS'*Strain*RotS;
  StrainF=RotD'*Strain1*RotD;
  
  RSxx(i)=StrainF(1,1);
  RSxy(i)=StrainF(1,2);
  RSxz(i)=StrainF(1,3);
  RSyy(i)=StrainF(2,2);
  RSyz(i)=StrainF(2,3);
  RSzz(i)=StrainF(3,3);
  
  %RSxx(i)=Strain(1,1);
  %RSxy(i)=Strain(1,2);
  %RSxz(i)=Strain(1,3);
  %RSyy(i)=Strain(2,2);
  %RSyz(i)=Strain(2,3);
  %RSzz(i)=Strain(3,3);
%  fprintf(Fid2,'%15.8f %15.8f %13.8f %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n',sx,sy,sz,RSxx,RSxy,RSxz,RSyy,RSyz,RSzz);
end
figure
min_s=-3*10e-6;
max_s= 3*10e-6;
index_max=RSyz>max_s;
index_min=RSyz<min_s;
index_nan=isnan(RSyz);
color_index=jet(128);
RRR=fix((RSyz-min_s)./((max_s-min_s)./128)+1);
RRR(index_max)=128;
RRR(index_min)=1;
RRR(index_nan)=1;
RR(:,1)=color_index(RRR,1);
RR(:,2)=color_index(RRR,2);
RR(:,3)=color_index(RRR,3);
scatter(sx,sy,20,RR,'filled')
%figure
%scatter(sx,sy,RSyz,'filled')
%fclose(Fid2);
%figure(15); quiver3(sx,sy,sz,Usumx,Usumy,Usumz,'r')
%figure(16); quiver(sx,sy,Usumx,Usumy,'r'); grid on;
end
