function least_pole
Fid1=fopen('/home/sasajima/DATA/PH-IA-lpdF.dat','r');
%Fid2=fopen('/home/sasajima/DATA/OUTPUT/PH-IH-leastpole.dat','w');
dF=textscan(Fid1,'%f %f %f %f %f %f %f %f %f');
dF=cell2mat(dF);
fclose(Fid1);

p0=7958.68798828;
B=[dF(:,3),dF(:,4),dF(:,5)];
A=[dF(:,1),dF(:,2)];
I=eye(3);
%P=[dF(:,9)/p0];
%P=diag(P);
%iP=inv(P);
%M=B*iP*B';
whos
M=B*I*B';
iM=inv(M);
w=[dF(:,6)];
K=A'*iM*A;
iK=inv(K);
whos
dlt=(-iK)*A'*iM*w
v=(I*B'*iM*A*iK*A'*iM*(-I)*B'*iM)*w;
whos
end
