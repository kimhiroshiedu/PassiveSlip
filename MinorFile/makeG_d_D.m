function[n]=makeG_d_D(trixyz3,xy)

%change the value%
pr=0.25;
ss= 0;
ds= 1;
ts= 0;
%==================
 
  n=length(trixyz3)
  p=length(xy);

  sx(:,1)=zeros(p,1);
  sy(:,1)=xy(:,1);
  sz(:,1)=xy(:,2);
  %sUxyz=zeros(n,n,3);
  dS6=zeros(p,n,6); 
  dUxyz=zeros(p,n,3);

  rootname4='/home_tmp/sasajima/DATA/MAT/longtime/1S6/S6';
  rootname3='/home_tmp/sasajima/DATA/MAT/longtime/1Uxyz/Uxyz';
  extension='.dat';
    
 for i=1:n%numbers of fault%
  
 di=i
 w=num2str(i);

  x=[trixyz3(i,1,1),trixyz3(i,1,2),trixyz3(i,1,3)];
  y=[trixyz3(i,2,1),trixyz3(i,2,2),trixyz3(i,2,3)];
  z=[trixyz3(i,3,1),trixyz3(i,3,2),trixyz3(i,3,3)];
  
  [U] = CalcTriDisps(sx, sy, sz, x, y, z, pr, ss, ts, ds);
  dUxyz(:,i,1)=U.x;
  dUxyz(:,i,2)=U.y;
  dUxyz(:,i,3)=-U.z;
  
  [S] = CalcTriStrains(sx, sy, sz, x, y, z, pr, ss, ts, ds);
  dS6(:,i,1)=S.xx;
  dS6(:,i,2)=S.xy;
  dS6(:,i,3)=S.xz;
  dS6(:,i,4)=S.yy;
  dS6(:,i,5)=S.yz;
  dS6(:,i,6)=S.zz;

  %dSss(:,i)=c1^2*Sxx+2*c1*s1*Sxy+s1^2*Syy;
  %dSsd(:,i)=c2*(c1*s1*(Syy-Sxx)+(c1^2-s1^2)*Sxy)+s2*(c1*Sxz+s1*Syz);
  %dSsn(:,i)=-s2(i).*(c1(i).*s1(i).*(Syy(:,i)-Sxx(:,i))+((c1(i)).^2-(s1(i)).^2).*Sxy(:,i))+c2(i).*(c1(i).*Sxz(:,i)+s1(i).*Syz(:,i));
  %dSdd(:,i)=c2^2*(s1^2*Sxx-2*s1*c1*Sxy+c1^2*Syy)+2*c2*s2*(c1*Syz-s1*Sxz)+s2^2*Szz;
  %dSdn(:,i)=c2(i).*s2(i).*(Szz(:,i)-(s1(i)).^2.*Sxx(:,i)+2.*s1(i).*c1(i).*Sxy(:,i)-(c1(i)).^2.*Syy(:,i))+((c2(i)).^2-(s2(i)).^2).*(c1(i).*Syz(:,i)-s1(i).*Sxz(:,i));
  %dSnn(:,i)=s2^2*(s1^2*Sxx-2*s1*c1*Sxy+c1^2*Syy)-2*c2*s2*(c1*Syz-s1*Sxz)+c2^2*Szz;

  filename4= [rootname4,w,extension];
  filename3= [rootname3,w,extension];

  S6i(:,1,1:6)=dS6(:,i,1:6);
  Uxyzi(:,1,1:3)=dUxyz(:,i,1:3);

  save(filename4,'S6i','-v7.3');
  save(filename3,'Uxyzi','-v7.3');
end
