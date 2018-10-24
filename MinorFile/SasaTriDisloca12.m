function SasaTriDisloca12
 %[triC,tri3,xy,tri]=Sasa_make_trixyz02;
 %save('/home/sasajima/DATA/MAT/triC_2','triC');
 %save('/home/sasajima/DATA/MAT/tri3_2','tri3');
 %save('/home/sasajima/DATA/MAT/tri_2','tri');
 %save('/home/sasajima/DATA/MAT/xy_2','xy');
 load('/home/sasajima/DATA/MAT/triC_2','triC');
 load('/home/sasajima/DATA/MAT/tri3_2','tri3');
 %[DistC]=distance_xyz(triC);
 %save('/home/sasajima/DATA/MAT/DistC','DistC');
 load('/home/sasajima/DATA/MAT/DistC','DistC');
 %[sSsn,sSdn,sUxyz]=makeG_s_03(triC,tri3);
 %save('/home/sasajima/DATA/MAT/sSsn_2','sSsn');
 %save('/home/sasajima/DATA/MAT/sSdn_2','sSdn');
 %save('/home/sasajima/DATA/MAT/sUxyz_2','sUxyz');
 %[dSsn,dSdn,dUxyz]=makeG_d_03(triC,tri3);
 %save('/home/sasajima/DATA/MAT/dSsn_2','dSsn');
 %save('/home/sasajima/DATA/MAT/dSdn_2','dSdn');
 %save('/home/sasajima/DATA/MAT/dUxyz_2','dUxyz');
 load('/home/sasajima/DATA/MAT/sSsn_2','sSsn');
 load('/home/sasajima/DATA/MAT/sSdn_2','sSdn');
 load('/home/sasajima/DATA/MAT/dSsn_2','dSsn');
 load('/home/sasajima/DATA/MAT/dSdn_2','dSdn');
 %[Slip,define]=defineSlip1(triC);
 %save('/home/sasajima/DATA/MAT/Slip_2','Slip');
 %save('/home/sasajima/DATA/MAT/define_2','define');
 load('/home/sasajima/DATA/MAT/Slip_2','Slip');
 load('/home/sasajima/DATA/MAT/define_2','define');
 %[A_Ssn,A_Sdn]=I_CalcStrain(Slip,sSsn,sSdn,dSsn,dSdn);
 %save('/home/sasajima/DATA/MAT/A_Ssn_2','A_Ssn');
 %save('/home/sasajima/DATA/MAT/A_Sdn_2','A_Sdn');
 %load('/home/sasajima/DATA/MAT/A_Ssn_2','A_Ssn');
 %load('/home/sasajima/DATA/MAT/A_Sdn_2','A_Sdn');
 %figure(11); quiver(triC(:,1),-triC(:,2),A_Ssn,A_Sdn,5,'r')

 %[R_Slip]=I_CalcSlip(A_Ssn,A_Sdn,sSsn,dSdn,Slip);
 %save('/home/sasajima/DATA/MAT/R_Slip_2','R_Slip');
 %load('/home/sasajima/DATA/MAT/I_Slip_2','I_Slip');
 %figure(12); quiver(triC(:,1),-triC(:,2),R_Slip(:,1),R_Slip(:,2),5,'g')
 
 n=length(Slip);
 %%[I_Slip]=make_I_Slip(Slip);
 %%save('/home/sasajima/DATA/MAT/I_Slip','I_Slip')
 %%load('/home/sasajima/DATA/MAT/I_Slip','I_Slip')
 R_Slip=Slip;
 SUM_Ssn=zeros(n,1);
 SUM_Sdn=zeros(n,1);
 SUM_Slip=Slip;
 
 for i=1:1                                                                       
 %figure(105); quiver(triC(:,1),triC(:,2),Slip(:,1),Slip(:,2),5,'b')
 %figure(105); quiver(triC(:,1),-triC(:,2),diag(sSsn),diag(dSdn),5,'b')
 %drawnow
 [R_Ssn,R_Sdn,SUM_Ssn,SUM_Sdn]=I2_CalcStrain12(R_Slip,SUM_Ssn,SUM_Sdn,sSsn,sSdn,dSsn,dSdn);
 %figure(101); quiver(triC(:,1),-triC(:,2),R_Ssn,R_Sdn,5,'r')
 %hold on
 %triplot(tri,xy(:,1),-xy(:,2))
 drawnow

 [R_Slip,K_Slip,C_Strain,SUM_Slip,sumsSsn,sumdSdn]=I2_CalcSlip12(DistC,R_Slip,SUM_Slip,R_Ssn,R_Sdn,sSsn,sSdn,dSsn,dSdn,Slip,define); 
 figure(1); quiver(triC(:,1),-triC(:,2),R_Slip(:,1),R_Slip(:,2),5,'g')
 figure(2); quiver(triC(:,1),-triC(:,2),K_Slip(:,1),K_Slip(:,2),5,'g')
 figure(3); quiver(triC(:,1),-triC(:,2),C_Strain(:,1),C_Strain(:,2),5,'r')
 S=sumsSsn';
 D=sumdSdn';
 %figure(4); quiver(triC(:,1),-triC(:,2),S,D,'r')
 %hold on
 %triplot(tri,xy(:,1),-xy(:,2))
 drawnow
 i
 pause(3)
 end
 
 %figure(102); quiver(triC(:,1),-triC(:,2),SUM_Ssn,SUM_Sdn,5,'r')
 %figure(2); quiver(triC(:,1),-triC(:,2),SUM_Slip(:,1),SUM_Slip(:,2),5,'g')
 
 save('/home/sasajima/DATA/MAT/IF_Ssn_12','SUM_Ssn');
 save('/home/sasajima/DATA/MAT/IF_Sdn_12','SUM_Sdn');
 %load('/home/sasajima/DATA/MAT/IF_Ssn_12','I_Ssn');
 %load('/home/sasajima/DATA/MAT/IF_Sdn_12','I_Sdn');
 save('/home/sasajima/DATA/MAT/IF_Slip_12','SUM_Slip');
 %%save('/home/sasajima/DATA/MAT/Sum_Slip','Sum_Slip');
 %load('/home/sasajima/DATA/MAT/IF_Slip_2','I_Slip');
 %%load('/home/sasajima/DATA/MAT/Sum_Slip','Sum_Slip');
end
