function [sSsn,sSdn,dSsn,dSdn]=loatMAT(trixyzC);

  rootname1='/home_tmp/sasajima/DATA/MAT/boso/sSsn/sSsn';
  rootname2='/home_tmp/sasajima/DATA/MAT/boso/sSdn/sSdn';
  rootname3='/home_tmp/sasajima/DATA/MAT/boso/dSsn/dSsn';
  rootname4='/home_tmp/sasajima/DATA/MAT/boso/dSdn/dSdn';

  extension='.dat';

  n=length(trixyzC);
  
 for i=1:n;

  loatMAT=i

  w=num2str(i);

  filename1= [rootname1,w,extension];
  filename2= [rootname2,w,extension];
  filename3= [rootname3,w,extension];
  filename4= [rootname4,w,extension];

  load(filename1,'sSsni','-mat');
  load(filename2,'sSdni','-mat');
  load(filename3,'dSsni','-mat');
  load(filename4,'dSdni','-mat');

 sSsn(:,i)=sSsni(:,1);
 sSdn(:,i)=sSdni(:,1);
 dSsn(:,i)=dSsni(:,1);
 dSdn(:,i)=dSdni(:,1);

 end
end
