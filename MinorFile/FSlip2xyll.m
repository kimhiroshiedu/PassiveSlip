function [xyFSlip]=FSlip2xyll(F_Slip,sitaS)
xyFSlip(:,1)=F_Slip(:,1).*cos(sitaS(:))+F_Slip(:,2).*sin(sitaS(:));
xyFSlip(:,2)=-F_Slip(:,1).*sin(sitaS(:))+F_Slip(:,2).*cos(sitaS(:));
end
