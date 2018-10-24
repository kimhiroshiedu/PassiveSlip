function calcUxyz[xyz,sUxyz,dUxyz,Slip]
%by Ryohei Sasajima
%last 2014/02/02
m=size(xyz,1);
n=size(Slip,1);

G(1:m,1:n)=sUxyz(1:m,1:n,1);
G(1:m,n+1:2*n)=dUxyz(1:m,1:n,1);

G(m+1:2*m,1:n)=sUxyz(1:m,1:n,2);
G(m+1:2*m,n+1:2*n)=sUxyz(1:m,1:n,2);

G(2*m+1:3*m,1:n)=sUxyz(1:m,1:n,3);
G(2*m+1:3*m,n+1:2*n)=sUxyz(1:m,1:n,3);

SSlip(1:n,1)=Slip(1:n,1);
SSlip(n+1:2*n,1)=Slip(1:n,2);

Uxyz(1:3*m,1)=G.*SSlip;

Fid1=fopen
Fid2=fopen

 for i=1:m
   fprintf(Fid1,'%12.6f %12.6f %12.6f %12.6f',xyz(i,1),xyz(i,2),Uxyz(i,1),Uxyz(m+i,1));
   fprintf(Fid2,'%12.6f %12.6f',xyz(i,1),xyz(i,2),Uxyz(2*m+i,1));
 end

end
