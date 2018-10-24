function RotNormal
%countour-clockwise rotation matrix
Z=[cos(sitaz),-sin(sitaz),0;sin(sitaz),cos(sitaz),0;0,0,1];
Y=[cos(sitay),0,sin(sitay);0,1,0;-sin(sitay),0,cos(sitay)];
X=[1,0,0;0,cos(sitax),-sin(sitax);0,sin(sitax),cos(sitax)];
end
