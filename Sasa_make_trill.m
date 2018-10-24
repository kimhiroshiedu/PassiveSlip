%====================================================
function [triC,tri3,tri,sll]=Sasa_make_trill
%====================================================
Fid=fopen('/home/sasajima/Dropbox/yellow/Kanto_bound2.txt','r');
bound=textscan(Fid,'%f %f');
fclose(Fid);
bound=cell2mat(bound);
%====================================================
Fid=fopen('/home/sasajima/Dropbox/yellow/PHdepth20130906.txt','r');
dep_main=textscan(Fid,'%f %f %f');
fclose(Fid);
dep_main=cell2mat(dep_main);
%====================================================
Fid=fopen('/home/sasajima/Dropbox/yellow/sg_collect.txt','r');
dep_sub=textscan(Fid,'%f %f %f');
fclose(Fid);
dep_sub=cell2mat(dep_sub);
%===================================================
%E=TriScatteredInterp(dep_main(:,1),dep_main(:,2),dep_main(:,3),'natural');
%F=TriScatteredInterp(dep_sub(:,1),dep_sub(:,2),dep_sub(:,3),'natural');
%min_lon=min(bound(:,1)); max_lon=max(bound(:,1));
%min_lat=min(bound(:,2)); max_lat=max(bound(:,2));
%figure
%plot(bound(:,1),bound(:,2),'r')
%hold on
%n=0;
%int_mesh=3000;
%while n<int_mesh
%  slat=(max_lat-min_lat).*rand(1)+min_lat;
%  slon=(max_lon-min_lon).*rand(1)+min_lon;
%  ID=inpolygon(slon,slat,bound(:,1),bound(:,2));
%  if ID==1
%    n=n+1;
%    s.lat(n)=slat;
%    s.lon(n)=slon;
    %s.dep(n)=F(slon,slat)-E(slon,slat);
    %if rem(n,round(int_mesh/10))==1;
    %  plot3(s.lon,s.lat,s.dep,'.')
    %  pause(.1)
    %end
%  end
%end
%====================================================
Fid=fopen('/home/sasajima/Dropbox/yellow/Kanto_mesh2.dat','r');
Tll=textscan(Fid,'%f %f %f');
fclose(Fid);
Tll=cell2mat(Tll);
s.lon(:,1)=Tll(:,1);
s.lat(:,1)=Tll(:,2);
sll=Tll;
%ll=[s.lon,s.lat];
%plot(s.lon,s.lat,'.')
%====================================================
tri = delaunay(s.lon,s.lat);
%====================================================
ntri=length(tri);
%nn=0;
%Fid=fopen('/home/sasajima/DATA/XYtriC.dat','w');
%Fid2=fopen('/home/sasajima/DATA/XYtri3.dat','w');
E=TriScatteredInterp(dep_main(:,1),dep_main(:,2),dep_main(:,3),'natural');
F=TriScatteredInterp(dep_sub(:,1),dep_sub(:,2),dep_sub(:,3),'natural');



%tri3=zeros(ntri,3,3);
%triC=zeros(ntri,3);
nn=0;
 for n=1:ntri
  
  slon1=s.lon(tri(n,1));
  slat1=s.lat(tri(n,1));
  sdep1=F(slon1,slat1)+E(slon1,slat1);
  
  slon2=s.lon(tri(n,2));
  slat2=s.lat(tri(n,2));
  sdep2=F(slon2,slat2)+E(slon2,slat2);
  
  slon3=s.lon(tri(n,3));
  slat3=s.lat(tri(n,3));
  sdep3=F(slon3,slat3)+E(slon3,slat3);
 
  if sdep1<0.1
    sdep1=0;
  else
  end
  
  if sdep2<0.1
    sdep2=0;
  else
  end

  if sdep3<0.1
    sdep3=0;
  else
  end

  glon=(slon1+slon2+slon3)/3;
  glat=(slat1+slat2+slat3)/3;
  gdep=(sdep1+sdep2+sdep3)/3;

  IN=inpolygon(glon,glat,bound(:,1),bound(:,2));

  if gdep==0
  elseif IN==0
  else
  nn=nn+1;
   triC(nn,:)=[glon,glat,gdep];
   tri3(nn,1,1)=slon1;
   tri3(nn,1,2)=slon2;
   tri3(nn,1,3)=slon3;
   tri3(nn,2,1)=slat1;
   tri3(nn,2,2)=slat2;
   tri3(nn,2,3)=slat3;
   tri3(nn,3,1)=sdep1;
   tri3(nn,3,2)=sdep2;
   tri3(nn,3,3)=sdep3;
  end
 end

  %fprintf(Fid,'%5d %11.5f %11.5f %11.5f\n',n,gx,gy,gdep);
  %fprintf(Fid2,'%5d %11.5f %11.5f %11.5f\n',n,sx1,sy1,sz1,sx2,sy2,sz2,sx3,sy3,sz3);
  %if ID==1
  %  nn=nn+1;
  %  s.tri(nn,:)=tri(n,:);
  %end
  %whos

%whos
%figure(47);
%plot(bound(:,1),bound(:,2),'r')
%hold on
%fclose(Fid);
%fclose(Fid2);
end
