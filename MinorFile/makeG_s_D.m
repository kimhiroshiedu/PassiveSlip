function[n]=makeG_s_D(trixyz3,xy)

%change the value%
pr=0.25;
ss= 1;
ds= 0;
ts= 0;
%==================
 
  n=length(trixyz3)
  p=length(xy);

  sx(:,1)=zeros(p,1);
  sy(:,1)=xy(:,1);
  sz(:,1)=xy(:,2);
  %sUxyz=zeros(n,n,3);
  sS6=zeros(p,n,6); 
  sUxyz=zeros(p,n,3);

  rootname4='/home_tmp/sasajima/DATA/MAT/longtime/2usS/sS';
  rootname3='/home_tmp/sasajima/DATA/MAT/longtime/2usUxyz/sUxyz';
  extension='.dat';
    
 for i=1:n%numbers of fault%
  
 si=i
 w=num2str(i);

  x=[trixyz3(i,1,1),trixyz3(i,1,2),trixyz3(i,1,3)];
  y=[trixyz3(i,2,1),trixyz3(i,2,2),trixyz3(i,2,3)];
  z=[trixyz3(i,3,1),trixyz3(i,3,2),trixyz3(i,3,3)];
  
  [U] = CalcTriDisps(sx, sy, sz, x, y, z, pr, ss, ts, ds);
  sUxyz(:,i,1)=U.x;
  sUxyz(:,i,2)=U.y;
  sUxyz(:,i,3)=-U.z;
  
  [S] = CalcTriStrains(sx, sy, sz, x, y, z, pr, ss, ts, ds);
  sS6(:,i,1)=S.xx;
  sS6(:,i,2)=S.xy;
  sS6(:,i,3)=S.xz;
  sS6(:,i,4)=S.yy;
  sS6(:,i,5)=S.yz;
  sS6(:,i,6)=S.zz;

  %sSss(:,i)=c1^2*Sxx+2*c1*s1*Sxy+s1^2*Syy;
  %sSsd(:,i)=c2*(c1*s1*(Syy-Sxx)+(c1^2-s1^2)*Sxy)+s2*(c1*Sxz+s1*Syz);
  %sSsn(:,i)=-s2(i).*(c1(i).*s1(i).*(Syy(:,i)-Sxx(:,i))+((c1(i)).^2-(s1(i)).^2).*Sxy(:,i))+c2(i).*(c1(i).*Sxz(:,i)+s1(i).*Syz(:,i));
  %sSdd(:,i)=c2^2*(s1^2*Sxx-2*s1*c1*Sxy+c1^2*Syy)+2*c2*s2*(c1*Syz-s1*Sxz)+s2^2*Szz;
  %sSdn(:,i)=c2(i).*s2(i).*(Szz(:,i)-(s1(i)).^2.*Sxx(:,i)+2.*s1(i).*c1(i).*Sxy(:,i)-(c1(i)).^2.*Syy(:,i))+((c2(i)).^2-(s2(i)).^2).*(c1(i).*Syz(:,i)-s1(i).*Sxz(:,i));
  %sSnn(:,i)=s2^2*(s1^2*Sxx-2*s1*c1*Sxy+c1^2*Syy)-2*c2*s2*(c1*Syz-s1*Sxz)+c2^2*Szz;

  filename4= [rootname4,w,extension];
  filename3= [rootname3,w,extension];

  sS6i(:,1,1:6)=sS6(:,i,1:6);
  sUxyzi(:,1,1:3)=sUxyz(:,i,1:3);

  save(filename4,'sS6i','-v7.3');
  save(filename3,'sUxyzi','-v7.3');
end
