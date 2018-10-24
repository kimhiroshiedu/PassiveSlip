%====================================================
function [triC,tri3]=Sasa_make_trixyz
%====================================================
%Fid=fopen('/home/sasajima/DATA/XYbound.txt','r');
%bound=textscan(Fid,'%f %f');
%fclose(Fid);
%bound=cell2mat(bound);
%====================================================
px=-270000;
n=0;
for k=1:26
 px=px+20000;
 py=-20000;
 for j=1:16
  n=n+1;
  py=py+20000;
  sx(n)=px;
  sy(n)=py;
 end
end
%====================================================
tri = delaunay(sx,sy);
%====================================================
ntri=length(tri);
%nn=0;
%Fid=fopen('/home/sasajima/DATA/XYtriC.dat','w');
%Fid2=fopen('/home/sasajima/DATA/XYtri3.dat','w');
for n=1:ntri
  gx=mean(sx(tri(n,:)));  
  gy=mean(sy(tri(n,:)));
  gz=-gy*0.3;
  %ID=inpolygon(gx,gy,bound(:,1),bound(:,2));
  
  sx1=sx(tri(n,1));
  sy1=sy(tri(n,1));
  sz1=-sy1*0.3;
  
  sx2=sx(tri(n,2));
  sy2=sy(tri(n,2));
  sz2=-sy2*0.3;
  
  sx3=sx(tri(n,3));
  sy3=sy(tri(n,3));
  sz3=-sy3*0.3;
  
  tri3(n,:,:)=[sx1,sx2,sx3;sy1,sy2,sy3;sz1,sz2,sz3];
  triC(n,:)=[gx,gy,gz];
  
  %fprintf(Fid,'%5d %11.5f %11.5f %11.5f\n',n,gx,gy,gdep);
  %fprintf(Fid2,'%5d %11.5f %11.5f %11.5f\n',n,sx1,sy1,sz1,sx2,sy2,sz2,sx3,sy3,sz3);
  %if ID==1
  %  nn=nn+1;
  %  s.tri(nn,:)=tri(n,:);
  %end  
end
%whos
%figure
%plot(bound(:,1),bound(:,2),'r')
%hold on
%triplot(s.tri,sx,sy);
%fclose(Fid);
%fclose(Fid2);
end
