function SasaTriDislocaBari
 %[triC,tri3,xy,tri]=Sasa_make_trixyz02;
 %save('/home/sasajima/DATA/MAT/triC_2','triC','tri3');
 load('/home/sasajima/DATA/MAT/tri_23','triC','tri3');

 %[sSsn,sSdn,sUxyz]=makeG_s_03(triC,tri3);
 %save('/home/sasajima/DATA/MAT/sSsn_2','sSsn');
 %save('/home/sasajima/DATA/MAT/sSdn_2','sSdn');
 %save('/home/sasajima/DATA/MAT/sUxyz_2','sUxyz');
 %[dSsn,dSdn,dUxyz]=makeG_d_03(triC,tri3);
 %save('/home/sasajima/DATA/MAT/dSsn_2','dSsn');
 %save('/home/sasajima/DATA/MAT/dSdn_2','dSdn');
 %save('/home/sasajima/DATA/MAT/dUxyz_2','dUxyz');
 load('/home/sasajima/DATA/MAT/sSsn_23','sSsn');
 load('/home/sasajima/DATA/MAT/sSdn_23','sSdn');
 load('/home/sasajima/DATA/MAT/dSsn_23','dSsn');
 load('/home/sasajima/DATA/MAT/dSdn_23','dSdn');
 
 [Area]=AreaCalc(tri3);
 save('/home/sasajima/DATA/MAT/Area23','Area');
 %load('/home/sasajima/DATA/MAT/Area23','Area');
 %
 
 [Slip,define]=defineSlipAB(triC);
 %save('/home_tmp/sasajima/DATA/MAT/Slip_2','Slip');
 %save('/home_tmp/sasajima/DATA/MAT/define_2','define');
 %load('/home_tmp/sasajima/DATA/MAT/Slip_2','Slip');
 %load('/home_tmp/sasajima/DATA/MAT/define_2','define');
 figure(200); quiver(triC(:,1),-triC(:,2),Slip(:,1),-Slip(:,2),2,'g')
 %
 %MSlip=Slip;
 %[Mw,M0]=MwCalc(Area,MSlip);
 %Mw
 %M0
 
 [A1_Ssn,A1_Sdn]=I_CalcStrain(Slip,sSsn,sSdn,dSsn,dSdn);
 %save('/home_tmp/sasajima/DATA/MAT/A1_Ssn_2','A1_Ssn');
 %save('/home_tmp/sasajima/DATA/MAT/A1_Sdn_2','A1_Sdn');
 %load('/home_tmp/sasajima/DATA/MAT/A1_Ssn_2','A1_Ssn');
 %load('/home_tmp/sasajima/DATA/MAT/A1_Sdn_2','A1_Sdn');
 %figure(201); quiver(triC(:,1),triC(:,2),A1_Ssn,A1_Sdn,5,'b')
 
 [A1_Slip]=OutofAsperity01(A1_Ssn,A1_Sdn,Slip,sSsn,sSdn,dSsn,dSdn);
 %save('/home_tmp/sasajima/DATA/MAT/A1Slip','A1_Slip'); 
 %load('/home_tmp/sasajima/DATA/MAT/A1Slip','A1_Slip');
 figure(1); quiver(triC(:,1),-triC(:,2),A1_Slip(:,1),-A1_Slip(:,2),2,'g')
 MSlip=A1_Slip;
 [Mw,M0]=MwCalc(Area,MSlip);
 Mw
 M0
 

 
 [B1_Ssn,B1_Sdn,AAA1_Slip]=I_CalcStrainAB(A1_Slip,sSsn,sSdn,dSsn,dSdn);
 %save('/home_tmp/sasajima/DATA/MAT/B1Strain','B1_Ssn','B1_Sdn');
 %load('/home_tmp/sasajima/DATA/MAT/B1Strain','B1_Ssn','B1_Sdn');
 figure(602); quiver(triC(:,1),-triC(:,2),AAA1_Slip(:,1),-AAA1_Slip(:,2),2,'r')
 figure(202); quiver(triC(:,1),-triC(:,2),B1_Ssn,-B1_Sdn,5,'b')
 
 %
 [B1_Slip]=OutofAsperityAB(B1_Ssn,B1_Sdn,A1_Slip,sSsn,sSdn,dSsn,dSdn);
 %save('/home_tmp/sasajima/DATA/MAT/A1B1Slip','B1_Slip');
 %load('/home_tmp/sasajima/DATA/MAT/A1B1Slip','B1_Slip');
 figure(52); quiver(triC(:,1),-triC(:,2),B1_Slip(:,1)+AAA1_Slip(:,1),-B1_Slip(:,2)-AAA1_Slip(:,2),2,'r')
 %whos
 %
 A1B1_Slip(:,1:2)=A1_Slip(:,1:2)+B1_Slip(:,1:2);
 index_barri=(Slip(:,3)==2);
 A1B1_Slip(index_barri,1)=0;
 A1B1_Slip(index_barri,2)=0;
 figure(3); quiver(triC(:,1),-triC(:,2),A1B1_Slip(:,1),-A1B1_Slip(:,2),2,'g')
 MSlip=A1B1_Slip;
 %whos
 [Mw,M0]=MwCalc(Area,MSlip);
 Mw
 M0
 %{
 load('/home/sasajima/DATA/MAT/FSlip2','F_Slip');
 AA_Slip(:,1:2)=A1B1_Slip(:,1:2)+F_Slip(:,1:2);
 MSlip=F_Slip;
 %whos
 [Mw,M0]=MwCalc(Area,MSlip);
 Mw
 M0
 
 
 %figure(33); quiver(triC(:,1),triC(:,2),AA_Slip(:,1),AA_Slip(:,2),2,'g')
 %whos
 %afterSlip(:,1:2)=A1_Slip(:,1:2)-A1B1_Slip(:,1:2);
 %figure(503); quiver(triC(:,1),triC(:,2),afterSlip(:,1),afterSlip(:,2),2,'r')
 %}
 %F_Slip=AA_Slip;
 F_Slip=A1B1_Slip;
 %F_Slip=afterSlip;
 %MSlip=F_Slip;
 %[Mw,M0]=MwCalc(Area,MSlip);
 %Mw
 %M0
 
 [F_Ssn,F_Sdn]=I_CalcStrainF(F_Slip,sSsn,sSdn,dSsn,dSdn);
 %save('/home_tmp/sasajima/DATA/MAT/FStrain','F_Ssn','F_Sdn');
 %load('/home_tmp/sasajima/DATA/MAT/FStrain','F_Ssn','F_Sdn');
 figure(101); quiver(triC(:,1),-triC(:,2),F_Ssn,-F_Sdn,5,'b')
 
 %absO_Ssn=abs(O_Ssn);
 %absO_Sdn=abs(O_Sdn);
 %sumAOSsn=sum(O_Ssn,1)
 %sumAOSdn=sum(O_Sdn,1)

 %}
 Fid=fopen('./sASlipF.dat','w');
 %Fid2=fopen('./sAuBStrainF.dat','w');
 Fid3=fopen('./sAuBSlipF.dat','w');
 Fid4=fopen('./sAuBStrainF.dat','w');
 
 n=size(Slip,1);
 for i=1:n
 fprintf(Fid, '%11.3f %11.3f %9.5f %9.5f\n',triC(i,1), triC(i,2), A1_Slip(i,1), A1_Slip(i,2));
 %fprintf(Fid2, '%11.3f %11.3f %13.5e %13.5e\n',triC(i,1), triC(i,2), B1_Ssn(i), B1_Sdn(i));
 fprintf(Fid3, '%11.3f %11.3f %9.5f %9.5f\n',triC(i,1), triC(i,2), A1B1_Slip(i,1), A1B1_Slip(i,2));
 fprintf(Fid4, '%11.3f %11.3f %13.5e %13.5e\n',triC(i,1), triC(i,2), F_Ssn(i), F_Sdn(i));
 end
 
 fclose(Fid);
 fclose(Fid2);
 fclose(Fid3);
 %fclose(Fid4);
 %{
 [FO_Ssn,FO_Sdn]=OutofAsperity11(F_Ssn,F_Sdn,Slip);
 %save('/home/sasajima/DATA/MAT/FOStrain','FO_Ssn','FO_Sdn'); 
 %load('/home/sasajima/DATA/MAT/FOStrain','FO_Ssn','FO_Sdn');
 absFO_Ssn=abs(FO_Ssn);
 absFO_Sdn=abs(FO_Sdn);
 sumF0Ssn=sum(absFO_Ssn,1)
 sumF0Sdn=sum(absFO_Sdn,1)
 %}
end
