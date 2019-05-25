function rho = RhoHelp_fun(e_ground,thetas,nsv,na)
%Polarization
% 10.7.18
%=======================
er = e_ground;
A = dot(na,nsv);
if abs(A)==1   % V
    a0=er;
    a1=a0*sin(thetas);
    a2=sqrt(er-cos(thetas)^2);
    rho = (a1-a2)./(a1+a2);
elseif (A==0)   % H
    a0=1;
    a1=a0*sin(thetas);
    a2=sqrt(er-cos(thetas)^2);
    rho = (a1-a2)./(a1+a2);
else
    rho =0;
    disp('error')
end
