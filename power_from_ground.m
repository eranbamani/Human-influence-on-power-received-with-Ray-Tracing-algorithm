function Pr_ref = power_from_ground(pr,pt,nsm,e_ground,npolv,k_wavenum,L_LOS,Gd)

p_ref = [pr(1),pr(2),2*pr(3)];
% 
rIv = p_ref-pt; 
L_ref= norm(rIv);
% 
costetagi = dot(rIv,nsm(2,:))/(L_ref*norm(nsm(2,:)));
thetagi = acos(costetagi);
theta_G = abs(thetagi-pi/2);
% 
rho = RhoHelp_fun(e_ground,theta_G,nsm(2,:),npolv);
% 
Delta_pi = k_wavenum*(L_ref-L_LOS);
Pr_ref = (rho).*Gd.*exp(-1i*Delta_pi)./L_ref; %Received Reflected Power [w]

end

