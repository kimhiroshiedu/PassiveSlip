function[si]=makeG_s_Q(trixyzC,trixyz3,sitaS,sitaD,normVec)

%change the value%
pr=0.25;
ss= 1;
ds= 0;
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
         ssn='sSsn/sSsn';
         sdn='sSdn/sSdn';
  extension='.dat';
    
 for i=1:n%numbers of fault%

 si=i
 w=num2str(i);

  x=[trixyz3(i,1,1),trixyz3(i,1,2),trixyz3(i,1,3)];
  y=[trixyz3(i,2,1),trixyz3(i,2,2),trixyz3(i,2,3)];
  z=[trixyz3(i,3,1),trixyz3(i,3,2),trixyz3(i,3,3)];
  
  %[U] = CalcTriDisps(sx, sy, sz, x, y, z, pr, ss, ts, ds);
  %sUxyz(:,i,1)=U.x;
  %sUxyz(:,i,2)=U.y;
  %sUxyz(:,i,3)=-U.z;
  [S] = CalcTriStrains(sx, sy, sz, x, y, z, pr, ss, ts, ds);
  Sxx(:,1)=S.xx;
  Sxy(:,1)=S.xy;
  Sxz(:,1)=S.xz;
  Syy(:,1)=S.yy;
  Syz(:,1)=S.yz;
  Szz(:,1)=S.zz;

  %sSss(:,i)=c1^2*Sxx+2*c1*s1*Sxy+s1^2*Syy;
  %sSsd(:,i)=c2*(c1*s1*(Syy-Sxx)+(c1^2-s1^2)*Sxy)+s2*(c1*Sxz+s1*Syz);
  sSsn(:,1)=-s2(i).*(c1(i).*s1(i).*(Syy(:,1)-Sxx(:,1))+((c1(i)).^2-(s1(i)).^2).*Sxy(:,1))+c2(i).*(c1(i).*Sxz(:,1)+s1(i).*Syz(:,1));
  %sSdd(:,i)=c2^2*(s1^2*Sxx-2*s1*c1*Sxy+c1^2*Syy)+2*c2*s2*(c1*Syz-s1*Sxz)+s2^2*Szz;
  sSdn(:,1)=c2(i).*s2(i).*(Szz(:,1)-(s1(i)).^2.*Sxx(:,1)+2.*s1(i).*c1(i).*Sxy(:,1)-(c1(i)).^2.*Syy(:,1))+((c2(i)).^2-(s2(i)).^2).*(c1(i).*Syz(:,1)-s1(i).*Sxz(:,1));
  %sSnn(:,i)=s2^2*(s1^2*Sxx-2*s1*c1*Sxy+c1^2*Syy)-2*c2*s2*(c1*Syz-s1*Sxz)+c2^2*Szz;

  filename1= [rootname,ssn,w,extension];
  filename2= [rootname,sdn,w,extension]; 

  sSsni(:,1)=sSsn(:,1);
  sSdni(:,1)=sSdn(:,1);

  save(filename1,'sSsni','-v7.3');
  save(filename2,'sSdni','-v7.3'); 
 end 
end
