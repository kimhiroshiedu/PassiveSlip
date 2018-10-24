function SasaTriDisloca14
 %[triC,tri3]=Sasa_make_trixyz04;
 %save('/home/sasajima/DATA/MAT/tri_curve1','triC','tri3');
 load('/home/sasajima/DATA/MAT/tri_curve1','triC','tri3');
 
 trixyzC=triC;
 trixyz3=tri3;
 
 %[sitaS,sitaD,normVec]=strike_dip(trixyzC,trixyz3);
 %save('/home_tmp/sasajima/DATA/MAT/sita_curve1','sitaS','sitaD','normVec');
 load('/home_tmp/sasajima/DATA/MAT/sita_curve1','sitaS','sitaD','normVec');
 
 %triC(:,3)
 
 %sitaD
 %
 
 [sSsn,sSdn]=makeG_s_Q(trixyzC,trixyz3,sitaS,sitaD,normVec);
 save('/home/sasajima/DATA/MAT/sSsn_curve1','sSsn');
 save('/home/sasajima/DATA/MAT/sSdn_curve1','sSdn');
 %save('/home/sasajima/DATA/MAT/sUxyz_18','sUxyz');
 [dSsn,dSdn]=makeG_d_Q(trixyzC,trixyz3,sitaS,sitaD,normVec);
 save('/home/sasajima/DATA/MAT/dSsn_curve1','dSsn');
 save('/home/sasajima/DATA/MAT/dSdn_curve1','dSdn');
 %save('/home/sasajima/DATA/MAT/sUxyz_18','sUxyz');
 %
 
 %{
 load('/home/sasajima/DATA/MAT/sSsn_curve1','sSsn');
 load('/home/sasajima/DATA/MAT/sSdn_curve1','sSdn');
 load('/home/sasajima/DATA/MAT/dSsn_curve1','dSsn');
 load('/home/sasajima/DATA/MAT/dSdn_curve1','dSdn');

 %figure(16);
 %quiver(trixyzC(:,1),trixyzC(:,2),normVec(:,1),normVec(:,2),1,'b')
 whos
 %triC

 [Slip,define]=defineSlip2(triC);
 %save('/home/sasajima/DATA/MAT/Slip_22','Slip');
 %save('/home/sasajima/DATA/MAT/define_22','define');
 %load('/home/sasajima/DATA/MAT/Slip_2','Slip');
 %load('/home/sasajima/DATA/MAT/define_2','define');
 figure(200); quiver(triC(:,1),-triC(:,2),Slip(:,1),-Slip(:,2),2,'g')

 [A_Ssn,A_Sdn]=I_CalcStrain(Slip,sSsn,sSdn,dSsn,dSdn);
 %save('/home/sasajima/DATA/MAT/A_Ssn_22','A_Ssn');
 %save('/home/sasajima/DATA/MAT/A_Sdn_22','A_Sdn');
 %load('/home/sasajima/DATA/MAT/A_Ssn_2','A_Ssn');
 %load('/home/sasajima/DATA/MAT/A_Sdn_2','A_Sdn');
 whos
 figure(201); quiver(triC(:,1),-triC(:,2),A_Ssn,-A_Sdn,5,'b')
 
 [F_Slip]=OutofAsperity01(A_Ssn,A_Sdn,Slip,sSsn,sSdn,dSsn,dSdn);
 %save('/home/sasajima/DATA/MAT/FSlip2','F_Slip'); 
 %load('/home/sasajima/DATA/MAT/FSlip','F_Slip');

 figure(1); quiver(triC(:,1),-triC(:,2),F_Slip(:,1),-F_Slip(:,2),2,'g')
 
 [F_Ssn,F_Sdn]=I_CalcStrainF(F_Slip,sSsn,sSdn,dSsn,dSdn);
 %save('/home/sasajima/DATA/MAT/FStrain2','F_Ssn','F_Sdn'); 
 %load('/home/sasajima/DATA/MAT/FStrain','F_Ssn','F_Sdn');
 figure(101); quiver(triC(:,1),-triC(:,2),F_Ssn,-F_Sdn,5,'b')
 
 whos
 
 %absO_Ssn=abs(O_Ssn);
 %absO_Sdn=abs(O_Sdn);
 %sumAOSsn=sum(O_Ssn,1)
 %sumAOSdn=sum(O_Sdn,1)
 
 [FO_Ssn,FO_Sdn]=OutofAsperity11(F_Ssn,F_Sdn,Slip);
 %save('/home/sasajima/DATA/MAT/FOStrain2','FO_Ssn','FO_Sdn'); 
 %load('/home/sasajima/DATA/MAT/FOStrain','FO_Ssn','FO_Sdn');
 
 %
 Fid=fopen('./2ACSlip.dat','w');
 Fid2=fopen('./2ACStrain1.dat','w');
 Fid3=fopen('./2ACSlipF.dat','w');
 Fid4=fopen('./2ACStrainF.dat','w');
 
 n=size(Slip,1);
 
 for i=1:n
 fprintf(Fid, '%11.3f %11.3f %9.5f %9.5f\n',triC(i,1), triC(i,2), Slip(i,1), Slip(i,2));
 fprintf(Fid2, '%11.3f %11.3f %13.5e %13.5e\n',triC(i,1), triC(i,2), A_Ssn(i), A_Sdn(i));
 fprintf(Fid3, '%11.3f %11.3f %9.5f %9.5f\n',triC(i,1), triC(i,2), F_Slip(i,1), F_Slip(i,2));
 fprintf(Fid4, '%11.3f %11.3f %13.5e %13.5e\n',triC(i,1), triC(i,2), F_Ssn(i), F_Sdn(i));
 end
 fclose(Fid);
 fclose(Fid2);
 fclose(Fid3);
 fclose(Fid4);
 %
 absFO_Ssn=abs(FO_Ssn);
 absFO_Sdn=abs(FO_Sdn);
 sumF0Ssn=sum(absFO_Ssn,1)
 sumF0Sdn=sum(absFO_Sdn,1)
 %}
end
