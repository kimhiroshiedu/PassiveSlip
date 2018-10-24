function [ll]=makell

W=138;
E=148;
dEW=0.1%(degree);
S=33;
N=48;
dSN=0.075%(degree);

nEW=(E-W)./dEW+1
nSN=(N-S)./dSN+1

m=nEW.*nSN

ll=zeros(m,2);

mm=0;

EW=W-dEW;

 for i=1:nEW;
 EW=EW+dEW;
 SN=S-dSN;
  for j=1:nSN; 
  SN=SN+dSN;
  mm=mm+1;
  ll(mm,1:2)=[EW,SN];
  end
 end

end
