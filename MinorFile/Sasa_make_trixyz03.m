%====================================================
function [triC,tri3]=Sasa_make_trixyz03
%====================================================
%Fid=fopen('/home/sasajima/DATA/XYbound.txt','r');
%bound=textscan(Fid,'%f %f');
%fclose(Fid);
%bound=cell2mat(bound);
%====================================================
px=-300000;
w=0;
for k=1:11
 px=px+50000;
 py=-50000;
 for j=1:7
  w=w+1;
  py=py+50000;
  sx(w)=px;
  sy(w)=py;
 end
end
%====================================================
tri = delaunay(sx,sy);
%====================================================
ntri=length(tri);
%nn=0;
%Fid=fopen('/home/sasajima/DATA/XYtriC.dat','w');
%Fid2=fopen('/home/sasajima/DATA/XYtri3.dat','w');
tri3=zeros(ntri,3,3);
triC=zeros(ntri,3);
for n=1:ntri
  %ID=inpolygon(gx,gy,bound(:,1),bound(:,2));
  
  sx1=sx(tri(n,1));
  sy1=sy(tri(n,1));
  sz1=sy1*0.8;
  
  sx2=sx(tri(n,2));
  sy2=sy(tri(n,2));
  sz2=sy2*0.8;
  
  sx3=sx(tri(n,3));
  sy3=sy(tri(n,3));
  sz3=sy3*0.8;
  
  gx=(sx1+sx2+sx3)/3;                
  gy=(sy1+sy2+sy3)/3;
  gz=gy*0.8;

  tri3(n,1,1)=sx1;
  tri3(n,1,2)=sx2;
  tri3(n,1,3)=sx3;
  tri3(n,2,1)=sy1;
  tri3(n,2,2)=sy2;
  tri3(n,2,3)=sy3;
  tri3(n,3,1)=sz1;
  tri3(n,3,2)=sz2;
  tri3(n,3,3)=sz3;
  triC(n,:)=[gx,gy,gz];
  %whos
  %fprintf(Fid,'%5d %11.5f %11.5f %11.5f\n',n,gx,gy,gdep);
  %fprintf(Fid2,'%5d %11.5f %11.5f %11.5f\n',n,sx1,sy1,sz1,sx2,sy2,sz2,sx3,sy3,sz3);

    
end
%whos
%figure
%plot(bound(:,1),bound(:,2),'r')
%hold on
%triplot(tri,sx,sy)
%fclose(Fid);
%fclose(Fid2);
end
