function[di]=makeG_d_Q(trixyzC,trixyz3,sitaS,sitaD,normVec)

%change the value%
pr=0.25;
ss= 0;
ds= 1;
ts= 0;
%==================
 
  n=length(trixyzC);

  sx=trixyzC(:,1)+normVec(:,1);
  sy=trixyzC(:,2)+normVec(:,2);
  sz=trixyzC(:,3)+normVec(:,3);
  c1=cos(sitaS);
  s1=sin(sitaS);
  c2=cos(sitaD);
  s2=sin(sitaD);

   rootname='/home_tmp/sasajima/DATA/GreenF/PACtest2PACtest/'
         dds='dSsn/dSsn';
         ddn='dSdn/dSdn';
  extension='.dat';
    
 for i=1:n%numbers of fault%
  
 di=i
 w=num2str(i);

  x=[trixyz3(i,1,1),trixyz3(i,1,2),trixyz3(i,1,3)];
  y=[trixyz3(i,2,1),trixyz3(i,2,2),trixyz3(i,2,3)];
  z=[trixyz3(i,3,1),trixyz3(i,3,2),trixyz3(i,3,3)];
  
  %[U] = CalcTriDisps(sx, sy, sz, x, y, z, pr, ss, ts, ds);
  %dUxyz(:,i,1)=U.x;
  %dUxyz(:,i,2)=U.y;
  %dUxyz(:,i,3)=-U.z;
  
  [S] = CalcTriStrains(sx, sy, sz, x, y, z, pr, ss, ts, ds);
  Sxx(:,1)=S.xx;
  Sxy(:,1)=S.xy;
  Sxz(:,1)=S.xz;
  Syy(:,1)=S.yy;
  Syz(:,1)=S.yz;
  Szz(:,1)=S.zz;

  %dSss(:,i)=c1^2*Sxx+2*c1*s1*Sxy+s1^2*Syy;
  %dSsd(:,i)=c2*(c1*s1*(Syy-Sxx)+(c1^2-s1^2)*Sxy)+s2*(c1*Sxz+s1*Syz);
  dSsn(:,i)=-s2(i).*(c1(i).*s1(i).*(Syy(:,1)-Sxx(:,1))+((c1(i)).^2-(s1(i)).^2).*Sxy(:,1))+c2(i).*(c1(i).*Sxz(:,1)+s1(i).*Syz(:,1));
  %dSdd(:,i)=c2^2*(s1^2*Sxx-2*s1*c1*Sxy+c1^2*Syy)+2*c2*s2*(c1*Syz-s1*Sxz)+s2^2*Szz;
  dSdn(:,i)=c2(i).*s2(i).*(Szz(:,1)-(s1(i)).^2.*Sxx(:,1)+2.*s1(i).*c1(i).*Sxy(:,1)-(c1(i)).^2.*Syy(:,1))+((c2(i)).^2-(s2(i)).^2).*(c1(i).*Syz(:,1)-s1(i).*Sxz(:,1));
  %dSnn(:,i)=s2^2*(s1^2*Sxx-2*s1*c1*Sxy+c1^2*Syy)-2*c2*s2*(c1*Syz-s1*Sxz)+c2^2*Szz;

  filename1= [rootname,dds,w,extension];
  filename2= [rootname,ddn,w,extension]; 

  dSsni(:,1)=dSsn(:,1);
  dSdni(:,1)=dSdn(:,1);

  save(filename1,'dSsni','-v7.3');
  save(filename2,'dSdni','-v7.3'); 
 end 
end
