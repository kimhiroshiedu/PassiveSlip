function diviSV

 Fid3=fopen('/home/sasajima/Dropbox/yellow/addPAC_aft311.txt','r');
 tmeca3=textscan(Fid3,'%f %f %f %f %f %f %f %f %f %f %f %f %f');
 fclose(Fid3);
 meca3=cell2mat(tmeca3);

 Fid4=fopen('/home/sasajima/Dropbox/yellow/addPAC_before311.txt','r');
 tmeca4=textscan(Fid4,'%f %f %f %f %f %f %f %f %f %f %f %f %f');
 fclose(Fid4);
 meca4=cell2mat(tmeca4);

 Fid1=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEPAC_thrust97_11Mw9.txt','r');
 tmeca1=textscan(Fid1,'%f %f %f %f %f %f %f %f %f %f %f %f %f');
 fclose(Fid1);
 meca1=cell2mat(tmeca1);

 Fid2=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEPAC_thrust11Mw9_140215.txt','r');
 tmeca2=textscan(Fid2,'%f %f %f %f %f %f %f %f %f %f %f %f %f');
 fclose(Fid2);
 meca2=cell2mat(tmeca2);

 si1=size(meca1,1);
 si2=size(meca2,1);
 si3=size(meca3,1);
 si4=size(meca4,1);

 meca(1:si1,1:13)=meca1;
 meca((si+1):(si1+si2),1:13)=meca2;
 meca((si1+si2+1):(si1+si2+si3),1:13)=meca3;
 meca((si1+si2+si3+1):(si1+si2+si3+si4),1:13)=meca4;
 
 n=size(meca,1)

nj=27;

cumNM0=zeros(nj,1)
cumSV=zeros(nj,1)
NWcumSV=zeros(nj,1)
NWcount=zeros(nj,1)

   jj=0;

 for j=36:49:0.5

   jj=jj+1;

   for i=1:n

     if (meca(i,2)>j)&&(meca(i,2)<j+1)

      M0=meca(i,10)^(meca(i,11)-7);             
      cumNM0(jj,1)=cumNM0(jj,1)+M0;
      cumSV(jj,1)=cumSV(jj,1)+(meca(i,4)+90)*M0;
      NWcumSV(jj,1)=NWcumSV(jj,1)+meca(i,4)+90;
      NWcount(jj,1)NWconut(jj,1)+1;
     else
     end

   end

 avewSV(1:nj,1)=cumSV(1:nj,1)./cumNM0(1:nj,1);%clockwise from North, (toward trench-->);
 avenwSV(1:nj,1)=NWcumSV(1:nj,1)./NWcount(1:nj,1);%non weighted

 fprintf
 fprintf

 end

end

