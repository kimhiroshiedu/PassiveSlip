function SasaTriDisloca11
 [triC,tri3]=Sasa_make_trixyz03;
 save('/home/sasajima/DATA/MAT/triC_3','triC');
 save('/home/sasajima/DATA/MAT/tri3_3','tri3');
 %load('/home/sasajima/DATA/MAT/triC_3','triC');
 %load('/home/sasajima/DATA/MAT/tri3_3','tri3');
 [DistC]=distance_xyz(triC);
 save('/home/sasajima/DATA/MAT/DistC','DistC');
 %load('/home/sasajima/DATA/MAT/DistC','DistC');
 %[sSsn,sSdn,sUxyz]=makeG_s_04(triC,tri3);
 %save('/home/sasajima/DATA/MAT/sSsn_2','sSsn');
 %save('/home/sasajima/DATA/MAT/sSdn_2','sSdn');
 %save('/home/sasajima/DATA/MAT/sUxyz_2','sUxyz');
 %[dSsn,dSdn,dUxyz]=makeG_d_04(triC,tri3);
 %save('/home/sasajima/DATA/MAT/dSsn_2','dSsn');
 %save('/home/sasajima/DATA/MAT/dSdn_2','dSdn');
 %save('/home/sasajima/DATA/MAT/dUxyz_2','dUxyz');
 load('/home/sasajima/DATA/MAT/sSsn_2','sSsn');
 load('/home/sasajima/DATA/MAT/sSdn_2','sSdn');
 load('/home/sasajima/DATA/MAT/dSsn_2','dSsn');
 load('/home/sasajima/DATA/MAT/dSdn_2','dSdn');
 [Slip,define]=defineSlip1(triC);
 save('/home/sasajima/DATA/MAT/Slip_3','Slip');
 save('/home/sasajima/DATA/MAT/define_3','define');
 %load('/home/sasajima/DATA/MAT/Slip_3','Slip');
 %load('/home/sasajima/DATA/MAT/define_3','define');
 [A_Ssn,A_Sdn]=I_CalcStrain(Slip,sSsn,sSdn,dSsn,dSdn);
 save('/home/sasajima/DATA/MAT/A_Ssn_3','A_Ssn');
 save('/home/sasajima/DATA/MAT/A_Sdn_3','A_Sdn');
 %load('/home/sasajima/DATA/MAT/A_Ssn_3','A_Ssn');
 %load('/home/sasajima/DATA/MAT/A_Sdn_3','A_Sdn');
 %%[I_Slip]=I_CalcSlip(A_Ssn,A_Sdn,sSsn,sSdn,dSsn,dSdn,Slip,define);
 %%save('/home/sasajima/DATA/MAT/I_Slip_3','I_Slip');
 %%load('/home/sasajima/DATA/MAT/I_Slip_3','I_Slip');
 
 %%n=length(Slip);
 %%[I_Slip]=make_I_Slip(Slip);
 %%save('/home/sasajima/DATA/MAT/I_Slip','I_Slip')
 %%load('/home/sasajima/DATA/MAT/I_Slip','I_Slip')
 I_Slip=Slip;
 for i=1:50
 %figure(105); quiver(triC(:,1),triC(:,2),Slip(:,1),Slip(:,2),5,'b')
 %figure(105); quiver(triC(:,1),triC(:,2),diag(sSsn),diag(dSdn),5,'b')
 %drawnow
 [I_Ssn,I_Sdn]=I2_CalcStrain(I_Slip,sSsn,sSdn,dSsn,dSdn);
 figure(103); quiver(triC(:,1),triC(:,2),I_Ssn,I_Sdn,5,'r')
 drawnow

 [I_Slip]=I2_CalcSlip2(DistC,I_Slip,A_Ssn,A_Sdn,I_Ssn,I_Sdn,sSsn,dSdn,Slip,define); 
 figure(100); quiver(triC(:,1),triC(:,2),I_Slip(:,1),I_Slip(:,2),5,'g')
 drawnow

 pause(1)
 end
 save('/home/sasajima/DATA/MAT/IF_Ssn_3','I_Ssn');
 save('/home/sasajima/DATA/MAT/IF_Sdn_3','I_Sdn');
 %load('/home/sasajima/DATA/MAT/IF_Ssn_3','I_Ssn');
 %load('/home/sasajima/DATA/MAT/IF_Sdn_3','I_Sdn');
 save('/home/sasajima/DATA/MAT/IF_Slip_3','I_Slip');
 %%save('/home/sasajima/DATA/MAT/Sum_Slip','Sum_Slip');
 %load('/home/sasajima/DATA/MAT/IF_Slip_3','I_Slip');
 %%load('/home/sasajima/DATA/MAT/Sum_Slip','Sum_Slip');
end
