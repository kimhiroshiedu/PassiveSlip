function circleSV

 Fid3=fopen('/home/sasajima/Dropbox/yellow/addPAC_aft311.txt','r');
 tmeca3=textscan(Fid3,'%f %f %f %f %f %f %f %f %f %f %f %f %f');
 fclose(Fid3);
 meca3=cell2mat(tmeca3);

 Fid4=fopen('/home/sasajima/Dropbox/yellow/addPAC_before311.txt','r');
 tmeca4=textscan(Fid4,'%f %f %f %f %f %f %f %f %f %f %f %f %f');
 fclose(Fid4);
 meca4=cell2mat(tmeca4);

 Fid1=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEPAC_th_befoMw9.txt','r');
 tmeca1=textscan(Fid1,'%f %f %f %f %f %f %f %f %f %f %f %f %f');
 fclose(Fid1);
 meca1=cell2mat(tmeca1);

 Fid2=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEPAC_th_aftMw9.txt','r');
 tmeca2=textscan(Fid2,'%f %f %f %f %f %f %f %f %f %f %f %f %f');
 fclose(Fid2);
 meca2=cell2mat(tmeca2);

 Fid5=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEwSVall.txt','w');
 
 Fid6=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEnwSVall.txt','w');

 Fid10=fopen('/home_tmp/sasajima/DATA/SVmesh.dat','r');
 tmblack=textscan(Fid10,'%f %f');
 fclose(Fid10);
 mblack=cell2mat(tmblack);

 Fid11=fopen('/home/sasajima/Dropbox/yellow/PAC_OKH_thrust.txt','r');
 tbound=textscan(Fid11,'%f %f');
 fclose(Fid11);
 bound=cell2mat(tbound);

 rr=size(mblack,1);
 s=0;

 for r=1:rr;

 Red=inpolygon(mblack(r,1),mblack(r,2),bound(:,1),bound(:,2));
 
   if Red==0
   else
    s=s+1
    meshll(s,1)=mblack(r,1);
    meshll(s,2)=mblack(r,2);
   end

 end%for r=1:rr

 si1=size(meca1,1);
 si2=size(meca2,1);
 si3=size(meca3,1);
 si4=size(meca4,1);

 meca(1:si1,1:13)=meca1;
 meca((si1+1):(si1+si2),1:13)=meca2;
 meca((si1+si2+1):(si1+si2+si3),1:13)=meca3;
 meca((si1+si2+si3+1):(si1+si2+si3+si4),1:13)=meca4;

 n=size(meca,1)

 nj=size(meshll,1)
 
 ll(1:nj,1:2)=meshll(1:nj,1:2);
 [xy]=ll2xy(ll);
 meshXY=xy;
 
 clearvars ll xy;

 ll(1:n,1:2)=meca(1:n,1:2);
 [xy]=ll2xy(ll);
 mecaXY=xy;

 thsh=50;%[km]

 for jj=1:nj
 jj
 cumNM0=0;
 cumSV=0;
 NWcumSV=0;
 NWcount=0;

     for i=1:n    

       dist=sqrt((meshXY(jj,1)-mecaXY(i,1))^2+(meshXY(jj,2)-mecaXY(i,2))^2);

       if dist<thsh
          
         M0=meca(i,10)^(meca(i,11)-7);
         cumNM0(1,1)=cumNM0(1,1)+M0;
         cumSV(1,1)=cumSV(1,1)+(meca(i,4)+90)*M0;
         NWcumSV(1,1)=NWcumSV(1,1)+meca(i,4)+90;
	 NWcount(1,1)=NWcount(1,1)+1;

       else
       end

     end

 if cumNM0==0;
  cumNM0=100;
 else
 end

 if NWcount==0;
  NWcount=1;
 else
 end

 avewSV(1,1)=cumSV(1,1)/cumNM0(1,1);%clockwise from North, (toward trench-->);
 avenwSV(1,1)=NWcumSV(1,1)/NWcount(1,1);%non weighted

 fprintf(Fid5,'%11.3f %11.3f %11.3f\n',meshll(jj,1),meshll(jj,2),avewSV(1,1)); 
 fprintf(Fid6,'%11.3f %11.3f %11.3f\n',meshll(jj,1),meshll(jj,2),avenwSV(1,1)); 
 

 end

end

