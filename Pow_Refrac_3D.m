function Pow_H_Refraction = Pow_Refrac_3D(pr,pt,po,lamda,k_wavenum,Gd,Pt)
%Refraction function
%ERAN BAMANI
%10.02.19
%-------------------------------------------------------------------
[vec,h,d_t_o,d_o_r,Pow_H_KED] = PowKED_3D(pr,pt,po,lamda,Gd);
alpah = h *((d_t_o+d_o_r)/(d_t_o*d_o_r));












end

