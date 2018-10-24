function [sUxyz,dUxyz]=loadMAT2(xyz,trixyz3);

  rootname='/home_tmp/sasajima/DATA/GreenF/PAC2test/';
  sU='sU';
  dU='dU';
  extension='.dat';

  m=size(xyz,1);
  n=size(trixyz3,1);
  
 for i=1:n;

  loadMAT=i

  w=num2str(i);

  filename1= [rootname,sU,w,extension];
  filename2= [rootname,dU,w,extension];

  load(filename1,'sUxyzi','-mat');
  load(filename2,'dUxyzi','-mat');

 sUxyz(:,i,1:3)=sSsni(:,1:3);
 dUxyz(:,i,1:3)=dSdni(:,1);

 end
end
