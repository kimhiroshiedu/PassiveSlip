function  [xyz]=makexyz
%by Ryohei Sasajima 
%last 2014/02/02
%====================================================
%{
Fid1=fopen('/home_tmp/sasajima/DATA/GEODETIC_DATA/coordinates_F3/sitelocate.txt','r');
geonet=textscan(Fid1,'%f %f %f');
fclose(Fid1);
geonet=cell2mat(geonet);
%}

Fid2=fopen('/home_tmp/sasajima/DATA/xyz01.dat','r');
xyz01=textscan(Fid2,'%f %f');
fclose(Fid2);
xyz01=cell2mat(xyz01);

xyz=xyz01;

end
