function SasaTriDislocaBoso
 
 %[triC,tri3,tri,sll]=Sasa_make_trill;
 %save('/home_tmp/sasajima/DATA/MAT/tri_Boso','triC','tri3','tri','sll');
 load('/home_tmp/sasajima/DATA/MAT/tri_Boso','triC','tri3','tri','sll');
 
 
 %[trixyzC,trixyz3,sxyz]=Sasa_ll2xyz(triC,tri3,sll);
 %save('/home_tmp/sasajima/DATA/MAT/trixyz_Boso','trixyzC','trixyz3','sxyz');
 load('/home_tmp/sasajima/DATA/MAT/trixyz_Boso','trixyzC','trixyz3','sxyz');


 %[sitaS,sitaD,normVec]=strike_dip(trixyzC,trixyz3);
 %save('/home_tmp/sasajima/DATA/MAT/sitaBoso','sitaS','sitaD','normVec');
 load('/home_tmp/sasajima/DATA/MAT/sitaBoso','sitaS','sitaD','normVec');

 n=size(sitaS,1)
 
 %figure(31); quiver(trixyzC(:,1),-trixyzC(:,2),sitaS,sitaD,1,'b')
 %
 %[sSsn,sSdn]=makeG_s_Q(trixyzC,trixyz3,sitaS,sitaD,normVec);
 %save('/home_tmp/sasajima/DATA/MAT/sSsn_Boso','sSsn','-v7.3');
 %save('/home_tmp/sasajima/DATA/MAT/sSdn_Boso','sSdn','-v7.3');
 %save('/home_tmp/sasajima/DATA/MAT/sxyz_Q11','sUxyz');
 %[dSsn,dSdn]=makeG_d_Q(trixyzC,trixyz3,sitaS,sitaD,normVec);
 %save('/home_tmp/sasajima/DATA/MAT/dSsn_Boso','dSsn','-v7,3');
 %save('/home_tmp/sasajima/DATA/MAT/dSdn_Boso','dSdn','-v7.3');
 %save('/home_tmp/sasajima/DATA/MAT/dUxyz_Q11','dUxyz');

 %{
 load('/home_tmp/sasajima/DATA/MAT/sSsn_Boso','sSsn');
 load('/home_tmp/sasajima/DATA/MAT/sSdn_Boso','sSdn');
 load('/home_tmp/sasajima/DATA/MAT/dSsn_Boso','dSsn');
 load('/home_tmp/sasajima/DATA/MAT/dSdn_Boso','dSdn');
 %}

 %
 [sSsn,sSdn,dSsn,dSdn]=loatMAT(trixyzC); 

 finish_loatMAT=0
 zeroV=zeros(n,1);
 %tri3(1243,3,:)
 whos
 %figure(22); quiver(trixyzC(:,1),-trixyzC(:,2),dSdn(:,1155),zeroV(:,1),5,'b')
 
 %
 %triC(1155,:)
 %
 [llF]=defineEllipsll;

 Velo(1,1)=0.0290;%(m/y)
 Velo(2,1)=0.0267;
 Velo(3,1)=0.0279;
 Velo(4,1)=0.0263;
 Velo(5,1)=0.0265;
 Velo(6,1)=0.0255;
 Velo(7,1)=0;

 %{
 Vec(1,1:3)=[139.350,35.300,325.0];%epilon,epilat,PlateVec crockwise from North [digree]
 Vec(2,1:3)=[139.680,35.140,336.5];
 Vec(3,1:3)=[139.950,35.660,329.4];
 Vec(4,1:3)=[140.160,34.680,333.6];
 Vec(5,1:3)=[139.970,34.920,335.0];
 Vec(6,1:3)=[140.595,34.700,332.5];
 Vec(7,1:3)=[140.500,35.075,346.0];
 %}

 %triC(1:10,:)
 %llF(1,:,:)
 %
 %[SVxy]=SV_ll2xy_boso(Vec)
 SVxy(1,1)=215.0./180.*pi;
 SVxy(2,1)=215.0./180.*pi;
 SVxy(3,1)=215.0./180.*pi;
 SVxy(4,1)=215.0./180.*pi;
 SVxy(5,1)=215.0./180.*pi;
 SVxy(6,1)=215.0./180.*pi;
 SVxy(7,1)=215.0./180.*pi;

 %
 [Slip]=defineSlipBoso(triC,sitaS,llF,Velo,SVxy);
 %finish_defineSlip=0
 %
 save('/home_tmp/sasajima/DATA/MAT/bososlip/sdSlip_boso','Slip');
 %load('/home_tmp/sasajima/DATA/MAT/bososlip/sdSlip_boso','Slip');
 %Slip
 figure(1); quiver(trixyzC(:,1),-trixyzC(:,2),Slip(:,1),Slip(:,2),5,'b')
 %
 whos

 %
 [A_Ssn,A_Sdn]=I_CalcStrain(Slip,sSsn,sSdn,dSsn,dSdn);
 finish_1stStrain=0
 %save('/home_tmp/sasajima/DATA/MAT/bososlip/A_Ssn_boso','A_Ssn','A_Sdn');
 %load('/home_tmp/sasajima/DATA/MAT/bososlip/A_Sdn_boso','A_Ssn','A_Sdn');
 figure(201); quiver(trixyzC(:,1),-trixyzC(:,2),A_Ssn,A_Sdn,5,'b')
 
 %whos

 [F_Slip]=OutofAsperity01(A_Ssn,A_Sdn,Slip,sSsn,sSdn,dSsn,dSdn);
 finish_F_Slip=0
 %save('/home_tmp/sasajima/DATA/MAT/bososlip/FSlip_boso','F_Slip'); 
 %load('/home_tmp/sasajima/DATA/MAT/bososlip/FSlip_boso','F_Slip');

 [xyFSlip]=FSlip2xyll(F_Slip,sitaS);
 finish_xyFSlip=0
 save('/home_tmp/sasajima/DATA/MAT/bososlip/xyFSlip_boso','xyFSlip');
 %
 load('/home_tmp/sasajima/DATA/MAT/bososlip/xyFSlip_boso','xyFSlip');
 figure(4); quiver(trixyzC(:,1),-trixyzC(:,2),xyFSlip(:,1),xyFSlip(:,2),5,'b')
 
 %Fid=fopen('./QxySlip11.dat','w');

 %n=size(Slip,1);
 %

 %{
[FSlipll]=vecxy2ll(trixyzC,xyFSlip)
 finish_FSlipll=0
 save('/home_tmp/sasajima/DATA/MAT/bososlip/FSlipll_boso','FSlipll');
 %load('/home_tmp/sasajima/DATA/MAT/bososlip/FSlipll_boso','FSlipll');
 figure(2); quiver(triC(:,1),triC(:,2),FSlipll(:,1),FSlipll(:,2),5,'b')
 %}

 %

 %
 ll=triC;

 [vnorm,v,NARFA]=SasaEular2Velo(ll)

 Vec=NARFA;
 [SVxy]=SV_ll2xy_boso(Vec)

%{
Value=SVxy(:,1).*180./pi;
sx=trixyzC(:,1);
sy=trixyzC(:,2);

min_s=35
max_s=36
index_max=Value>max_s;
index_min=Value<min_s;
index_nan=isnan(Value);
color_index=jet(128);
RRR=fix((Value-min_s)./((max_s-min_s)./128)+1);
RRR(index_max)=128;
RRR(index_min)=1;
RRR(index_nan)=1;
RR(:,1)=color_index(RRR,1);
RR(:,2)=color_index(RRR,2);
RR(:,3)=color_index(RRR,3);
figure(19);
scatter(sx,-sy,20,RR,'filled') 
%}

 %
 [AnormaryVec]=AnormaryVecCalc(vnorm,SVxy,xyFSlip)

  Fid2=fopen('/home_tmp/sasajima/DATA/AnormaryVec2.txt','w');
  m=length(ll);
  for i=1:m
   fprintf(Fid2, '%12.5f %12.5f %12.5f\n',ll(i,1),ll(i,2),AnormaryVec(i,1));
  end
  fclose(Fid2);
 %

 %for i=1:n
 %fprintf(Fid, '%11.3f %11.3f %9.5f %9.5f\n',trixyzC(i,1), trixyzC(i,2), xyFSlip(i,1), xyFSlip(i,2));
 %end
 %fclose(Fid);
 
 %{
 figure(1); quiver(trixyzC(:,1),-trixyzC(:,2),xyFSlip(:,1),-xyFSlip(:,2),2,'g')
 %hold on
 %triplot(tri,1000.*sxyz(:,1),1000.*sxyz(:,2));
 
 [F_Ssn,F_Sdn]=I_CalcStrainF(F_Slip,sSsn,sSdn,dSsn,dSdn);
 %save('/home_tmp/sasajima/DATA/MAT/FStrain','F_Ssn','F_Sdn'); 
 %load('/home_tmp/sasajima/DATA/MAT/FStrain','F_Ssn','F_Sdn');
 figure(401); quiver(trixyzC(:,1),-trixyzC(:,2),F_Ssn,F_Sdn,5,'b')
 
 %whos
 
 %absO_Ssn=abs(O_Ssn);
 %absO_Sdn=abs(O_Sdn);
 %sumAOSsn=sum(absO_Ssn,1)
 %sumAOSdn=sum(absO_Sdn,1)
 
 [FO_Ssn,FO_Sdn]=OutofAsperity11(F_Ssn,F_Sdn,Slip);
 %save('/home_tmp/sasajima/DATA/MAT/FOStrain','FO_Ssn','FO_Sdn'); 
 %load('/home_tmp/sasajima/DATA/MAT/FOStrain','FO_Ssn','FO_Sdn');
 absFO_Ssn=abs(FO_Ssn);
 absFO_Sdn=abs(FO_Sdn);
 sumF0Ssn=sum(absFO_Ssn,1)
 sumF0Sdn=sum(absFO_Sdn,1)
 %}
end
