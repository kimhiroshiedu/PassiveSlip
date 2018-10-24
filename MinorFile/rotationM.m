function rotationM

RsitaSX90=0.2; %countour clockwise from x-axis +pi [rad]%

RotS=[cos(RsitaSX90),-sin(RsitaSX90),0;sin(RsitaSX90),cos(RsitaSX90),0;0,0,1];

Rdip360=2.4; %countour clockwise from earth surface%

RotD=[1,0,0;0,cos(Rdip360),-sin(Rdip360);0,sin(Rdip360),cos(Rdip360)];

T0=[3,8,0;8,4,5;0,5,9];
  %StrainF=RotD'*RotS'*Strain*RotS*RotD;
T1=RotS'*T0*RotS
T2=RotD'*T1*RotD

end