function PassiveSlip
%by Ryohei Sasajima 
%final 2013/12/23

 
 %[triC,tri3,tri,sll]=make_test_trill;
 %save('/home_tmp/sasajima/DATA/PassiveSlip/PAC_test/tri','triC','tri3','tri','sll');
 load('/home_tmp/sasajima/DATA/PassiveSlip/PAC_test/tri','triC','tri3','tri','sll');
 
 
 %[trixyzC,trixyz3,sxyz]=trill2trixyz(triC,tri3,sll);
 %save('/home_tmp/sasajima/DATA/PassiveSlip/PAC_test/trixyz','trixyzC','trixyz3','sxyz');
 load('/home_tmp/sasajima/DATA/PassiveSlip/PAC_test/trixyz','trixyzC','trixyz3','sxyz');

 %Gtrixyz3=trixyz3;


 %[sitaS,sitaD,normVec]=strike_dip(trixyzC,trixyz3);
 %save('/home_tmp/sasajima/DATA/PassiveSlip/PAC_test/sita','sitaS','sitaD','normVec');
 load('/home_tmp/sasajima/DATA/PassiveSlip/PAC_test/sita','sitaS','sitaD','normVec');

 %trixyz3=Gtrixyz3;

 n=size(sitaS,1)
 
 %figure(31); quiver(trixyzC(:,1),-trixyzC(:,2),sitaS,sitaD,1,'b')
  
 [si]=makeG_s_Q(trixyzC,trixyz3,sitaS,sitaD,normVec);

 [di]=makeG_d_Q(trixyzC,trixyz3,sitaS,sitaD,normVec);

 %{

 [xyz]=makexyz;
 save('/home_tmp/sasajima/DATA/PassiveSlip/PAC_test/xyz','xyz');
 %load('/home_tmp/sasajima/DATA/PassiveSlip/PAC_test/xyz','xyz');

 [si]=makeG_s_O(xyz,trixyz3);

 [di]=makeG_d_O(xyz,trixyz3);

 [sUxyz,dUxyz]=loadMAT2(xyz);

 [sSsn,sSdn,dSsn,dSdn]=loadMAT(trixyzC);

 finish_loadMAT=0
 
 zeroV=zeros(n,1);
 %tri3(1243,3,:)
 whos
 %figure(22); quiver(trixyzC(:,1),-trixyzC(:,2),dSdn(:,1155),zeroV(:,1),5,'b')
 
 %
 %triC(1155,:)
 %
 [Slip]=defineSlipQ(triC,sitaS);
 save('/home_tmp/sasajima/DATA/MAT/sdSlip_Q1','Slip');
 %load('/home_tmp/sasajima/DATA/MAT/sdSlip_Q1','Slip');
 
 [xySlip]=Slip2xyll(Slip,sitaS);
 save('/home_tmp/sasajima/DATA/MAT/xySlip_Q1','xySlip');
 %load('/home_tmp/sasajima/DATA/MAT/xySlip_Q1','xySlip');
 %figure(103); quiver(triC(:,1),triC(:,2),xySlip(:,1),-xySlip(:,2),'b')
 figure(101); quiver(trixyzC(:,1),-trixyzC(:,2),xySlip(:,1),-xySlip(:,2),'g')
 
 ll(:,1:2)=triC(:,1:2);
 %[vnorm,v,NARFA]=SasaEular2Velo(ll);

 []=relavec2vecVxy(NARFA)

 [xyBslip]=BackSlipxy(mcmcC,triC,trixyzC,sitaS,accuV,vecVxy)

 [A_Ssn,A_Sdn]=I_CalcStrain(Slip,sSsn,sSdn,dSsn,dSdn);
 save('/home_tmp/sasajima/DATA/MAT/A_Ssn_2','A_Ssn','A_Sdn');
 %load('/home_tmp/sasajima/DATA/MAT/A_Sdn_2','A_Ssn','A_Sdn');
 figure(201); quiver(trixyzC(:,1),-trixyzC(:,2),A_Ssn,A_Sdn,5,'b')
 
 %whos
 
 [F_Slip]=OutofAsperity01(A_Ssn,A_Sdn,Slip,sSsn,sSdn,dSsn,dSdn);
 save('/home_tmp/sasajima/DATA/MAT/FSlip','F_Slip'); 
 %load('/home_tmp/sasajima/DATA/MAT/FSlip','F_Slip');

 [xyFSlip]=FSlip2xyll(F_Slip,sitaS);
 save('/home_tmp/sasajima/DATA/MAT/xyFSlip_Q1','xyFSlip');
 %load('/home_tmp/sasajima/DATA/MAT/xyFSlip_Q1','xyFSlip');
 
 
 %Fid=fopen('./QxySlip11.dat','w');

 
 %n=size(Slip,1);
 
 %for i=1:n
 %fprintf(Fid, '%11.3f %11.3f %9.5f %9.5f\n',trixyzC(i,1), trixyzC(i,2), xyFSlip(i,1), xyFSlip(i,2));
 %end
 %fclose(Fid);
 
 %
 figure(1); quiver(trixyzC(:,1),-trixyzC(:,2),xyFSlip(:,1),-xyFSlip(:,2),2,'g')
 %hold on
 %triplot(tri,1000.*sxyz(:,1),1000.*sxyz(:,2));
 
 [F_Ssn,F_Sdn]=I_CalcStrainF(F_Slip,sSsn,sSdn,dSsn,dSdn);
 save('/home_tmp/sasajima/DATA/MAT/FStrain','F_Ssn','F_Sdn'); 
 %load('/home_tmp/sasajima/DATA/MAT/FStrain','F_Ssn','F_Sdn');
 figure(401); quiver(trixyzC(:,1),-trixyzC(:,2),F_Ssn,F_Sdn,5,'b')
 
 %whos
 
 %absO_Ssn=abs(O_Ssn);
 %absO_Sdn=abs(O_Sdn);
 %sumAOSsn=sum(absO_Ssn,1)
 %sumAOSdn=sum(absO_Sdn,1)
 
 [FO_Ssn,FO_Sdn]=OutofAsperity11(F_Ssn,F_Sdn,Slip);
 save('/home_tmp/sasajima/DATA/MAT/FOStrain','FO_Ssn','FO_Sdn'); 
 %load('/home_tmp/sasajima/DATA/MAT/FOStrain','FO_Ssn','FO_Sdn');
 absFO_Ssn=abs(FO_Ssn);
 absFO_Sdn=abs(FO_Sdn);
 sumF0Ssn=sum(absFO_Ssn,1)
 sumF0Sdn=sum(absFO_Sdn,1)
 %}
end
