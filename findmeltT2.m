function [T_m] = findmeltT2(P)
out = fliplr(SF_PhaseLines('Ih','water1'));
phasei = out;
clear out
out = fliplr(SF_PhaseLines('III','water1'));
phaseiii  = out;
clear out
out = fliplr(SF_PhaseLines('V','water1'));
phasev = out;
clear out
out = fliplr(SF_PhaseLines('VI','water1'));
phasevi = out;
compil = [phasei;phaseiii;phasev;phasevi];

[~,ind] = unique(compil(:,1));

for i = 1:length(P)
T_m(i,1)=interp1(compil(ind,1),compil(ind,2),P(i),"linear",'extrap');
end
end
