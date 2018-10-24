function [xySlip]=Slip2xyll(Slip,sitaS)
xySlip(:,1)=Slip(:,1).*cos(sitaS(:))+Slip(:,2).*sin(sitaS(:));
xySlip(:,2)=-Slip(:,1).*sin(sitaS(:))+Slip(:,2).*cos(sitaS(:));
end
