function [Pr_tot] = PowRec_fun(L_ref,L_LOS,k_wavenum,rho_T,Gd)
%Power Receive
% 11.7.18
%=======================
Delta_pi = k_wavenum*(L_ref-L_LOS);
Pr_tot = (rho_T).*Gd.*exp(-1i*Delta_pi)./L_ref; %Received Power from Reflected [w]
end
