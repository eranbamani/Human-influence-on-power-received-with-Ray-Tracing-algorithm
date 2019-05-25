function [vec,h,d_t_o,d_o_r,Pow_H_KED,L] = PowKED_3D(pr,pt,po,lamda,Gd)
%Diffraction function
%ERAN BAMANI
%19.11.18
%-------------------------------------------------------------------
%Data
%vector calculation
rv1= po-pt; %vector Tx to obj
rv2= pr-po; %vector obj to Rx
rv3= pr-pt; %vector Rx to Tx
rv4=pt-po;
d_t_o= norm(rv1); % length from Tx to obj
d_o_r= norm(rv2); % length from obj to Rx
%----------------------
%Distance 
up2=norm(cross(rv3,rv4));
down2=norm(rv3);
h=up2/down2;
t=-(dot(rv4,rv3))./norm(rv3)^2;
vec=pt+(rv3).*t;
%----------------------
%Diffraction constant L(v)
v=h*sqrt((2*d_t_o*d_o_r)./(lamda*(d_t_o+d_o_r)));  %fresnel zone
if (-0.8<=v)&&(v<0)
    L=0.5-0.62*v;
elseif (0<=v)&&(v<=1)
    L=0.5*exp(-0.95*v);
elseif (1<=v)&&(v<=2.4)
    L=0.4-sqrt(0.1184-(0.38-0.1*v)^2);
elseif (v>2.4)
    L=0.225/v;
end
%----------------------
phase = (pi*2*(d_t_o+d_o_r))/lamda;
%Pow Diff
Pow_H_KED= L*Gd*exp(-1j*phase);
L = 20.*log10(L);
end

