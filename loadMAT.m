function [sSsn,sSdn,dSsn,dSdn]=loadMAT(trixyzC);

  rootname='/home_tmp/sasajima/DATA/GreenF/PACTohoku2PACTohoku/';
  ss='sSsn/sSsn';
  sd='sSdn/sSdn';
  ds='dSsn/dSsn';
  dd='dSdn/dSdn';
  extension='.dat';

  n=length(trixyzC);
  
 for i=1:n;

  loatMAT=i

  w=num2str(i);

  filename1= [rootname,ss,w,extension];
  filename2= [rootname,sd,w,extension];
  filename3= [rootname,ds,w,extension];
  filename4= [rootname,dd,w,extension];

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
