function[si]=makeG_s_O(xyz,trixyz3)

%change the value%
pr=0.25;
ss= 1;
ds= 0;
ts= 0;
%==================
 
  m=size(xyz,1);
  n=size(trixyz3,1);

  sx=xyz(:,1);
  sy=xyz(:,2);
  sz=zeros(m,1);

   rootname='/home_tmp/sasajima/DATA/GreenF/PAC2test/'
         sU='sU';
  extension='.dat';
    
 for i=1:n%numbers of fault%
  
 si=i
 w=num2str(i);

  x=[trixyz3(i,1,1),trixyz3(i,1,2),trixyz3(i,1,3)];
  y=[trixyz3(i,2,1),trixyz3(i,2,2),trixyz3(i,2,3)];
  z=[trixyz3(i,3,1),trixyz3(i,3,2),trixyz3(i,3,3)];
  
  [U] = CalcTriDisps(sx, sy, sz, x, y, z, pr, ss, ts, ds);
  sUxyz(:,1)=U.x;
  sUxyz(:,2)=U.y;
  sUxyz(:,3)=-U.z;
  
  %[S] = CalcTriStrains(sx, sy, sz, x, y, z, pr, ss, ts, ds);
  %Sxx(:,i)=S.xx;
  %Sxy(:,i)=S.xy;
  %Sxz(:,i)=S.xz;
  %Syy(:,i)=S.yy;
  %Syz(:,i)=S.yz;
  %Szz(:,i)=S.zz;

  filename1= [rootname,sU,w,extension];

  sUxyzi(:,1:3)=sUxyz(:,1:3);

  save(filename1,'sUxyzi','-v7.3'); 
 end 
end
