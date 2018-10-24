%====================================================
function [triC,tri3,tri,sll]=cmtmake_trill
%by Ryohei Sasajima 2013/12/23
%====================================================
%====================================================
Fid0=fopen('/home/sasajima/Dropbox/yellow/PACdepth201312.txt','r');
dep_main=textscan(Fid0,'%f %f %f');
fclose(Fid0);
dep_main=cell2mat(dep_main);
%====================================================
%====================================================
Fid1=fopen('/home/sasajima/Dropbox/yellow/PAC_all.txt','r');
bound=textscan(Fid1,'%f %f');
fclose(Fid1);
bound=cell2mat(bound);


Fid10=fopen('/home_tmp/sasajima/DATA/green.dat','r');
mgreen=textscan(Fid10,'%f %f');
fclose(Fid10);
mgreen=cell2mat(mgreen);
gg=size(mgreen,1);

s=0;

 for g=1:gg;

 Green=inpolygon(mgreen(g,1),mgreen(g,2),bound(:,1),bound(:,2));

   if Green==0
   else
     s=s+1;
     slon(s,1)=mgreen(g,1);
     slat(s,1)=mgreen(g,2);
   end
 end

 sall=s

 sll(:,1:2)=[slon(:,1),slat(:,1)];

%====================================================
tri = delaunay(slon,slat);
%====================================================
ntri=length(tri);

E=TriScatteredInterp(dep_main(:,1),dep_main(:,2),dep_main(:,3),'natural');%depth of interplate -km

nn=0;
 for n=1:ntri
  
  lon1=slon(tri(n,1));
  lat1=slat(tri(n,1));
  dep1=-E(lon1,lat1);%+km
  
  lon2=slon(tri(n,2));
  lat2=slat(tri(n,2));
  dep2=-E(lon2,lat2);
  
  lon3=slon(tri(n,3));
  lat3=slat(tri(n,3));
  dep3=-E(lon3,lat3);

  glon=(lon1+lon2+lon3)/3;
  glat=(lat1+lat2+lat3)/3;
  gdep=(dep1+dep2+dep3)/3;

   if gdep==0
   else
    nn=nn+1;
    triC(nn,1:3)=[glon,glat,gdep];
    tri3(nn,1,1)=lon1;
    tri3(nn,1,2)=lon2;
    tri3(nn,1,3)=lon3;
    tri3(nn,2,1)=lat1;
    tri3(nn,2,2)=lat2;
    tri3(nn,2,3)=lat3;
    tri3(nn,3,1)=dep1;
    tri3(nn,3,2)=dep2;
    tri3(nn,3,3)=dep3;
   end
 
 end

end
