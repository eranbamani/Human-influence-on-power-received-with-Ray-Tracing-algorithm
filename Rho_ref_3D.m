function rho_3d = Rho_ref_3D(eo,er_H,tg_alpha,pp,pt)
%Reflection Coefficients for PowRed_3D
%ERAN BAMANI
% 23.12.18
%=======================
er = eo \ er_H;
if pp(2)<pt(2) || pp(2)>pt(2) % H
    a0=1;
    a1=a0*sin(tg_alpha);
    a2=sqrt(er-cos(tg_alpha)^2);
    rho_3d = (a1-a2)./(a1+a2);
else
    rho_3d =1;
end








