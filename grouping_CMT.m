function GroupingCMT
%for grouping CMT mechanism of earthquakes
%by Ryohei Sasajima 
%final 2013/12/23

 
 %[triC,tri3,tri,sll]=cmtmake_trill;
 %save('/home_tmp/sasajima/DATA/PassiveSlip/PAC_cmt/tri','triC','tri3','tri','sll');
 load('/home_tmp/sasajima/DATA/PassiveSlip/PAC_cmt/tri','triC','tri3','tri','sll');
 

 %[trixyzC,trixyz3]=cmttrill2trixyz(triC,tri3,sll);
 %save('/home_tmp/sasajima/DATA/PassiveSlip/PAC_cmt/trixyz','trixyzC','trixyz3');
 load('/home_tmp/sasajima/DATA/PassiveSlip/PAC_cmt/trixyz','trixyzC','trixyz3'); 

 %trixyzC(1:10,:)
 %trixyz3(1:10,:,:)
 %pause

 %[sitaS,sitaD,normVec]=cmtstrike_dip(trixyzC,trixyz3);
 %save('/home_tmp/sasajima/DATA/PassiveSlip/PAC_cmt/sita','sitaS','sitaD','normVec');
 load('/home_tmp/sasajima/DATA/PassiveSlip/PAC_cmt/sita','sitaS','sitaD','normVec');


 ll(:,1:2)=triC(:,1:2);

 %sitaS(1:10,1)
 %sitaD(20:30,1)


 [vnorm,v,NARFA]=SasaEular2Velo(ll);

 [n]=groupCMT(NARFA,sitaS,sitaD,triC);

end

