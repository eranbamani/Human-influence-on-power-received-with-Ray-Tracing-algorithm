function [d1,d2,pp,Pow_REF] = PowRef_3D(pr,pt,po,Ro,eo,er_H_R,k_wavenum,Gd)
%Reflaction function
%ERAN BAMANI
%19.11.18
%-------------------------------------------------------------------
%Data
pr_xy=pr([1,2]);
pt_xy=pt([1,2]);
po_xy=po([1,2]);
%----------------------
%vector calculation for Los
rv = pr-pt;
L_LOS = norm(rv);
%----------------------
%vector calculation for point p
rtc= pt_xy-po_xy; %vector obj to Tx
rrc= pr_xy-po_xy; %vector obj to Rx
%----------------------
%normal of impact point
R_tc= rtc./norm(rtc); %Normalizing the vector
R_rc= rrc./norm(rrc); %Normalizing the vector
n_p = ((R_tc+R_rc)./(norm(R_tc+R_rc))); % normal for impact point
n_p_3D = n_p;
n_p_3D(3)=0;
R_pc=((R_tc+R_rc)./(norm(R_tc+R_rc))).*Ro; %Size of the vector between p-c
%----------------------
%Point calculation
pp_xy=R_pc+po_xy;
pp=pp_xy;
%Z point
AA= sqrt((pp(1)-pt(1))^2+(pp(2)-pt(2))^2);
BB= sqrt((pr(1)-pp(1))^2+(pr(2)-pp(2))^2);
z2=((pr(3)-pt(3))/AA)*((AA*BB)/(AA+BB));
z1= pr(3)-pt(3)-z2;
pp(3)=pt(3)+z1;
%----------------------
%Point%Distance calculation
vr_in=pp-pt;
d1=norm(vr_in);
vr_out=pr-pp;
d2=norm(vr_out);
L_ref3D=d1+d2;
%----------------------
%Angle - tg_alpha
rrp= pr-pp;
costetari = dot(rrp,n_p_3D)./(norm(rrp)*norm(n_p_3D));
thetari = acos(costetari);
tg_alpha = abs(thetari-pi/2);
%----------------------
%Pow REF
rho_3d = Rho_ref_3D(eo,er_H_R,tg_alpha,pp,pt);
Delta_pi_3D = k_wavenum*(L_ref3D-L_LOS); 
Pow_REF = (rho_3d).*Gd.*exp(-1i*Delta_pi_3D)./L_ref3D;


end

