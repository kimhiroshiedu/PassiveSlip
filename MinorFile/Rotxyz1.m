function Rotxyz1
%change the value%
Fid=fopen('/home/sasajima/DATA/Output/NChiba.dat','w');
sitax=31.3/180*pi;  %SlipVector of oceanic plate side counter-clockwise from North% 
sitay=2*pi-35.9/180*pi;  %lat
sitaz=140.5/180*pi;  %lon

%countour-clockwise rotation matrix
Z=[cos(sitaz),-sin(sitaz),0;sin(sitaz),cos(sitaz),0;0,0,1];
Y=[cos(sitay),0,sin(sitay);0,1,0;-sin(sitay),0,cos(sitay)];
X=[1,0,0;0,cos(sitax),-sin(sitax);0,sin(sitax),cos(sitax)];

EER=6371000.0;
sita0=0;
 for i=1:10000
 sita0=sita0+2*pi/10000;
  x1=EER*cos(sita0);
  y1=EER*sin(sita0);
  xyz1(:,i)=[x1;y1;0];
 end
  xyz2=X*xyz1;
  xyz3=Y*xyz2;
  xyzF=Z*xyz3;
  xyzF=xyzF';

 for i=1:10000
  xF=xyzF(i,1);
  yF=xyzF(i,2);
  zF=xyzF(i,3);
  fprintf(Fid,'%15.7f %15.7f %15.7f\n',xF,yF,zF);
 end
whos 
end
