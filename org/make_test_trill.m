%====================================================
function [triC,tri3,tri,sll]=make_test_trill
%by Ryohei Sasajima 2013/12/23
%====================================================
%====================================================
Fid0=fopen('/home/sasajima/Dropbox/yellow/PACdepth201312.txt','r');
dep_main=textscan(Fid0,'%f %f %f');
fclose(Fid0);
dep_main=cell2mat(dep_main);
%====================================================
Fid00=fopen('/home/sasajima/Dropbox/yellow/depth0_1402.txt','r');
dep_sub=textscan(Fid00,'%f %f %f');
fclose(Fid00);
dep_sub=cell2mat(dep_sub);
%===================================================
%====================================================
Fid1=fopen('/home/sasajima/Dropbox/yellow/PAC_blue.txt','r');
bound=textscan(Fid1,'%f %f');
fclose(Fid1);
bound=cell2mat(bound);

Fid2=fopen('/home/sasajima/Dropbox/yellow/PAC_blue.txt','r');
PACblue=textscan(Fid2,'%f %f');
fclose(Fid2);
PACblue=cell2mat(PACblue);

Fid3=fopen('/home/sasajima/Dropbox/yellow/PAC_green.txt','r');
PACgreen=textscan(Fid3,'%f %f');
fclose(Fid3);
PACgreen=cell2mat(PACgreen);

Fid4=fopen('/home_tmp/sasajima/DATA/blue.dat','r');
mblue=textscan(Fid4,'%f %f');
fclose(Fid4);
mblue=cell2mat(mblue);
bb=size(mblue,1);

Fid5=fopen('/home_tmp/sasajima/DATA/green.dat','r');
mgreen=textscan(Fid5,'%f %f');
fclose(Fid5);
mgreen=cell2mat(mgreen);
gg=size(mgreen,1);

Fid7=fopen('/home/sasajima/Dropbox/yellow/PAC_red.txt','r');
PACred=textscan(Fid7,'%f %f');
fclose(Fid7);
PACred=cell2mat(PACred);

Fid8=fopen('/home/sasajima/Dropbox/yellow/PAC_black.txt','r');
PACblack=textscan(Fid8,'%f %f');
fclose(Fid8);
PACblack=cell2mat(PACblack);

Fid9=fopen('/home_tmp/sasajima/DATA/red.dat','r');
mred=textscan(Fid9,'%f %f');
fclose(Fid9);
mred=cell2mat(mred);
rr=size(mred,1);

Fid10=fopen('/home_tmp/sasajima/DATA/black.dat','r');
mblack=textscan(Fid10,'%f %f');
fclose(Fid10);
mblack=cell2mat(mblack);
blbl=size(mblack,1);

Fid11=fopen('/home_tmp/sasajima/DATA/m_all.dat','r');
mall=textscan(Fid11,'%f %f');
fclose(Fid11);
mall=cell2mat(mall);
alal=size(mall,1);


s=0;

 for bl=1:blbl;

 Blue=inpolygon(mblack(bl,1),mblack(bl,2),PACblue(:,1),PACblue(:,2));

  if Blue==0
  else
   s=s+1;
   slon(s,1)=mblack(bl,1);
   slat(s,1)=mblack(bl,2);
  end
 end

 s

%{

 for b=1:bb;

 Blue=inpolygon(mblue(b,1),mblue(b,2),PACblue(:,1),PACblue(:,2));
 BRed=inpolygon(mblue(b,1),mblue(b,2),PACred(:,1),PACred(:,2));

  if (Blue==1)&&(BRed==0)
   s=s+1;
   slon(s,1)=mblue(b,1);
   slat(s,1)=mblue(b,2);
  else
  end
 end

 sblue=s-sred

 for g=1:gg;

  Green=inpolygon(mgreen(g,1),mgreen(g,2),PACgreen(:,1),PACgreen(:,2));
 GBlue=inpolygon(mgreen(g,1),mgreen(g,2),PACblue(:,1),PACblue(:,2));

  if (Green==1)&&(GBlue==0)
   s=s+1;
   slon(s,1)=mgreen(g,1);
   slat(s,1)=mgreen(g,2);
  else
  end
 end

 sgreen=s-sred-sblue

 for bl=1:blbl;

  Black=inpolygon(mblack(bl,1),mblack(bl,2),PACblack(:,1),PACblack(:,2));
 BlGreen=inpolygon(mblack(bl,1),mnlack(bl,2),PACgreen(:,1),PACgreen(:,2));

  if (Black==1)&&(BlGreen==0)
   s=s+1;
   slon(s,1)=mblack(bl,1);
   slat(s,1)=mblack(bl,2);
  else
  end
 end

 sblack=s-sred-sblue-sgreen


 for al=1:alal;

  All=inpolygon(mall(al,1),mall(al,2),bound(:,1),bound(:,2));
 AlBlack=inpolygon(mall(al,1),mall(al,2),PACblack(:,1),PACblack(:,2));

  if (All==1)&&(AlBlack==0)
   s=s+1;
   slon(s,1)=mall(al,1);
   slat(s,1)=mall(al,2);
  else
  end
 end

 soutblack=s-sred-sblue-sgreen-sblack
 sall=s
 pause

%}

 sll(:,1:2)=[slon(:,1),slat(:,1)];

%ll=[s.lon,s.lat];
%plot(s.lon,s.lat,'.')
%====================================================
tri = delaunay(slon,slat);
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
  
  lon1=slon(tri(n,1));
  lat1=slat(tri(n,1));
  dep1=-F(lon1,lat1)-E(lon1,lat1);
  
  lon2=slon(tri(n,2));
  lat2=slat(tri(n,2));
  dep2=-F(lon2,lat2)-E(lon2,lat2);
  
  lon3=slon(tri(n,3));
  lat3=slat(tri(n,3));
  dep3=-F(lon3,lat3)-E(lon3,lat3);
  %{
   if dep1<0.1
     dep1=0;
   else
   end
  
   if dep2<0.1
     dep2=0;
   else
   end

   if dep3<0.1
     dep3=0;
   else
   end
  %}
  glon=(lon1+lon2+lon3)/3;
  glat=(lat1+lat2+lat3)/3;
  gdep=(dep1+dep2+dep3)/3;

  IN=inpolygon(glon,glat,bound(:,1),bound(:,2));

   if gdep==0
   elseif IN==0
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

   Fid6=fopen('/home_tmp/sasajima/DATA/PAC_tri3.txt','w');

 for w=1:nn;
  for u=1:3;
   fprintf(Fid6,'%10.4f %9.4f %10.4f\n',tri3(w,1,u),tri3(w,2,u),tri3(w,3,u));
  end
 end
  
  fclose(Fid6);

  %fprintf(Fid7,'%5d %11.5f %11.5f %11.5f\n',n,gx,gy,gdep);
  %fprintf(Fid8,'%5d %11.5f %11.5f %11.5f\n',n,sx1,sy1,sz1,sx2,sy2,sz2,sx3,sy3,sz3);
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
