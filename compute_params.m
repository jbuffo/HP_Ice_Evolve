function [out,K] = compute_params(PT,phase_s)
Phase = phase_s;
if Phase==0;
    out=SeaFreeze(PT,'water1');
    K=.60;
elseif Phase==1;
    out=SeaFreeze(PT,'Ih');
    K= 2.4;
 
elseif Phase==2;
    out=SeaFreeze(PT,'II');
    K= 1.8;
 
elseif Phase==3;
    out=SeaFreeze(PT,'III');
    K= 1.1;
    
elseif Phase==5;
    out=SeaFreeze(PT,'V');
    K=1.5;
 
else %Phase==6;
    out=SeaFreeze(PT,'VI');
    K= 1.9;
 
end
end