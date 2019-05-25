% Asymmetrical
% Obj-Detect
% Eran BAMANI
% 21.11.18
%==========================================================================
clear all
format long
% Data
% Room size
x_Room= 100; %[m] x #half length of the room
y_Room=100; %[m] y #half width of the room
z_Room=100; %[m] z #half hight of the room
%----------------------
% Reflactions
icase_r=0; % 0-nmr, 1- nr,mr
%
mr =3; % Horizontal reflection
nr =3; % Vertical reflection
%
nmr=1; % total reflection: icase_r=0
%----------------------
% Reciver Data
xr = 0;
wr = 1.5; %[m]
hr = 1.06; %[m]
%----------------------
% Transmitter Data
xt = 3; %[m]
wt=1.5;  %[m]
ht=1.06;  %[m]
%----------------------
% Object Data
xoo = 1.5; %[m]
wo=0;  %[m]
ho=1.67;  %[m]
Ro=0.5;
%----------------------
% Impact surface
nsm(1,:) = [0,0,-1];  % Celing
nsm(2,:) = [0,0,1];   % Ground
nsm(3,:) = [0,-1,0];  % Wall right
nsm(4,:) = [0,1,0];   % Wall left
%----------------------
% Antenna Properties
Gt_db=1; %  Tx Gain [dB]
Gr_db=1; %  Rx Gain [dB]
Pt_dbm=50;   %  Power Tx [dBm]
f = 500; % [MHz] Freq.
npolv = [0,0,1]; %antenna polarization 
%----------------------
% Material Properties   % er = e_ / e0;
eo = 8.85e-12;
e_ground = 52; %Permitivity Ground
% e_celing=; %Permitivity Celing
% e_walls = ; %Permitivity Walls
% er_H_T = (eo*63.3-4.4e-10*1i)+0.7*(eo*80.3-1.5915e-11*1i); %Permitivity Object
er_H_R = 63.3 + 0.7*80.3;
% erV = [er1,er2,er3,er3]; %Permitivity Vector
%================================================
% Units
c=3E8; % speed of light [m/s]
%================================================
fr=f*1E6; % [Hz] 
Gt=10^(Gt_db/10); % [-]
Gr=10^(Gr_db/10); % [-]
Pt=(10^(Pt_dbm/10))/1e3;   %  Power Tx [W]
Gd=sqrt(Gt*Gr);
%
lamda = c/fr; % Wavelength in free space [m]
k_wavenum=(2*pi)/lamda; % Wavenumber 

%----------------------
% Geomertry
% Points Rx
yr = wr;%-y_Room; %From the wall to Center of the room
zr=  hr;%-z_Room; %From the Ground to Center of the room
pr=[xr,yr,zr];
% Points Tx
yt = wt;%-y_Room; %From the wall to Center of the room
zt = ht;%-z_Room; %From the Ground to Center of the room
pt=[xt,yt,zt]; %Transmitter Point
% Points Object
xo=xoo;%-x_Room;
zo=ho;%-z_Room;
% yo=wo-y_Room;
n_y = 121;
yo = linspace(0,3,n_y);
%----------------------
% Load Data
[data] = xlsread('Res_9_0.xlsx','D14:E4');
Pev=data(:,2)';
xev=data(:,1)';
%==========================================================================
for i = 1:length(yo)
    po=[xo,yo(i),zo]; %Object Point
    dis = abs(po(1)-pr(1)); %distance between xo to LOS
    %Power received
    PH = Pt*(lamda/(4*pi))^2; %%Antenna Aperture Constant / Free space Loss
%     if dis > Ro
        %---------------------
        % Radiation Model
        rv = pr-pt;
        L_LOS = norm(rv);
        P_LOS = Gd./L_LOS;
        Pr_ref = power_from_ground(pr,pt,nsm,e_ground,npolv,k_wavenum,L_LOS,Gd);
        Pr_radiation = P_LOS+Pr_ref;
%--------------influence of human body on power received---------------------
    % Reflaction
        [d1,d2,pp,Pow_H_REF] = PowRef_3D(pr,pt,po,Ro,eo,er_H_R,k_wavenum,Gd);
    % %----------------------
    % Knife Edge
        [vec,h,d_t_o,d_o_r,Pow_H_KED,L] = PowKED_3D(pr,pt,po,lamda,Gd);
    % %----------------------
    % Scattering
    %
    % %----------------------
    % Total Power Recive 
%     ####1---------------------------------
        Pr_total = PH*abs(Pr_radiation)^2;
        Pow_H_REF1 = PH*abs(Pow_H_REF)^2;
        Pow_H_KED1 = PH*abs(Pow_H_KED)^2;
%         Pow_H_REF1 = abs(Pow_H_REF)^2;
%         Pow_H_KED1 = abs(Pow_H_KED)^2;
%         Pow_H_REF1 = abs(Pow_H_REF);
%         Pow_H_KED1 = abs(Pow_H_KED);
        Pr_total_dBm(i) = 10.*log10(Pr_total+Pow_H_REF1+Pow_H_KED1);
%     ####2---------------------------------
        Pr_total2 = PH*abs(Pr_radiation+Pow_H_REF+Pow_H_KED)^2;
        Pr_total_dBm2(i) = 10.*log10(Pr_total2);
        
%     else
%         Pow_H_Refraction = Pow_Refrac_3D();
        %         [vec,h,d_t_o,d_o_r,Pow_H_KED] = PowKED_3D(pr,pt,po,lamda,k_wavenum,Gd,Pt);
%         Pow_H_KED = abs(Pow_H_KED);
%     end
end
%==========================================================================
% Total Power Recive in dBm
% Pr_tot_dBm = 10.*log10(Pr_Reflaction_TOT);  %Total Received Power  [dBm]
% Pow_REF_TOT_dBm = 10.*log10(Pow_REF_TOT);
% 

%==========================================================================
% ###Plot###
nFont=12;
nLine=2;
% 
figure(1)
plot(yo,Pr_total_dBm);
axis([0 3 -120 0])
hold on
% figure(2)
plot(yo,Pr_total_dBm2);
axis([0 3 -120 0])
hold on
% figure(3)
plot(xev,Pev)
axis([0 3 -120 0])
hold on



% xlabel('Distance [m]');
% ylabel('P_r [dBm]');
% title('Radiation&Reflaction');
% legend('Genral Ray','LOS')
% text(0.4,0.4,['Pr_total_dBm=',num2str(Pr_total_dBm)],'Units','normalized','Fontsize',nFont)
% text(0.8,0.6,['b[m]=',num2str(z_Room)],'Units','normalized','Fontsize',nFont)
% text(0.8,0.5,['mr=',num2str(mr)],'Units','normalized','Fontsize',nFont)
% text(0.8,0.4,['nr=',num2str(nr)],'Units','normalized','Fontsize',nFont)
% hold off
% set(gca,'Fontsize',nFont)
% grid on
% %========================================================================
% figure(2)
% plot(yo,P_LOS1);%,xr,Pr_LOS_dBm,'--');
% xlabel('Distance [m]');
% ylabel('P_r [dBm]');
% title('Genral Ray Tracing vs Meas');
% legend('Genral Ray Tracing','Experiment')
% hold on
% set(gca,'Fontsize',nFont)
% grid on
%==========================================================================

% plot
% yt_max=mr*2*y_Room+y_Room;  %for Axis
% zt_max=nr*2*z_Room+z_Room;  %for Axis
% y_start=-40;
% y_end=0;

%==========================================================================
% figure(3)
% plot(yt,zt,'X',yr,zr,'*',y_img,z_img,'o')
% grid on
% axis([-yt_max,yt_max,-zt_max,zt_max])
% hold on
% yy=[-z_Room,z_Room,z_Room,-z_Room,-z_Room];
% zz=[-y_Room,-y_Room,y_Room,y_Room,-y_Room];
% %
% yy1=[(2*j-1)*z_Room,(2*j+1)*z_Room,(2*j+1)*z_Room,(2*j-1)*z_Room,(2*j-1)*z_Room];
% zz1=[(2*i-1)*y_Room,(2*i-1)*y_Room,(2*i+1)*y_Room,(2*i+1)*y_Room,(2*i-1)*y_Room];
% %
% plot(yy,zz,yy1,zz1)
% hold off
%==========================================================================
% figure(4)
% text(0.1,0.4,['nrC=',num2str(nr_C)],'Units','normalized','Fontsize',nFont)
% text(0.1,0.3,['nrG=',num2str(nr_G)],'Units','normalized','Fontsize',nFont)
% text(0.1,0.2,['mrR=',num2str(mr_R)],'Units','normalized','Fontsize',nFont)
% text(0.1,0.1,['mrL=',num2str(mr_L)],'Units','normalized','Fontsize',nFont)
% 
% text(0.1,0.9,['Pr_L_O_S [dBm]=',num2str(Pr_LOS_dBm)],'Units','normalized','Fontsize',nFont)
% text(0.1,0.8,['Pr_t_o_t [dBm]=' ,num2str(Pr_tot_dBm)],'Units','normalized','Fontsize',nFont)

%==========================================================================

% figure(5)
% text(0.1,0.4,['Pr-LOS(dBm)=',num2str(Pr_LOS_dBm)],'Units','normalized','Fontsize',nFont)
% text(0.1,0.3,['Pr-tot(dBm)=',num2str(Pr_tot_dBm)],'Units','normalized','Fontsize',nFont)
% text(0.1,0.1,['Pr-LOS=',num2str(Pr_LOS)],'Units','normalized','Fontsize',nFont)
% text(0.1,0.2,['Pr-tot=',num2str(Pr_tot)],'Units','normalized','Fontsize',nFont)
%==========================================================================
% figure(6)
% plot(k_v,Pr_tot,'r',k_v,Pr_LOS,'--')
% grid on
% hold on
%==========================================================================
% Points for room
P2xv=[x_Room,x_Room,x_Room,x_Room,x_Room];
P1xv=-P2xv;
%
P1yv=[-y_Room,y_Room,y_Room,-y_Room,-y_Room];
P1zv=[-z_Room,-z_Room,z_Room,z_Room,-z_Room];
%
Q1xv=[-x_Room,x_Room];
Q1yv=[y_Room,y_Room];
Q1zv=[z_Room,z_Room];
% %----------------------
% figure(7) %Indoor & obj & Ref Ray
% plot3(P1xv,P1yv,P1zv,'k','LineWidth',3) %Front rectangle
% xlabel('x [m]')
% ylabel('y [m]')
% zlabel('z [m]')
% hold on
% plot3( P2xv, P1yv, P1zv,'k','LineWidth',3)%Rear rectangle
% plot3(Q1xv,Q1yv, Q1zv,'k','LineWidth',3) %Left up
% plot3(Q1xv,-Q1yv, Q1zv,'k','LineWidth',3)%Right up
% plot3(Q1xv,-Q1yv,-Q1zv,'k','LineWidth',3)%Right down
% plot3(Q1xv,Q1yv,-Q1zv,'k','LineWidth',3) %Left down
% % Tx&Rx
% plot3( pt(1),pt(2),pt(3),'O')
% plot3( [pt(1) pt(1)],[pt(2) pt(2)],[-z_Room pt(3)],'r','LineWidth',2)
% plot3(pr(1) ,pr(2),pr(3),'*')
% plot3( [pr(1) pr(1)],[pr(2) pr(2)],[-z_Room pr(3)],'r','LineWidth',2)
% % LOS
% plot3( [pt(1) pr(1)],[pt(2) pr(2)],[pt(3) pr(3)],'--','LineWidth',2)
% Reflection
%1
% plot3( 0,y_Room,0,'X')
% plot3( [pt(1) 0],[pt(2) y_Room],[pt(3) 0],'-.','LineWidth',0.8)
% plot3( [0 pr(1)],[y_Room pr(2)],[0 pr(3)],'-.','LineWidth',0.8)
% %2
% plot3( 0,-y_Room,0,'X')
% plot3( [pt(1) 0],[pt(2) -y_Room],[pt(3) 0],'-.','LineWidth',0.8)
% plot3( [0 pr(1)],[-y_Room pr(2)],[0 pr(3)],'-.','LineWidth',0.8)
% %3
% plot3( 0,0,z_Room,'X')
% plot3( [pt(1) 0],[pt(2) 0],[pt(3) z_Room],'-.','LineWidth',0.8)
% plot3( [0 pr(1)],[0 pr(2)],[z_Room pr(3)],'-.','LineWidth',0.8)
% %4
% plot3( 0,0,-z_Room,'X')
% plot3( [pt(1) 0],[pt(2) 0],[pt(3) -z_Room],'-.','LineWidth',0.8)
% plot3( [0 pr(1)],[0 pr(2)],[-z_Room pr(3)],'-.','LineWidth',0.8)
% grid on
% Obj
%-------------------------------------------
% plot3(xo,yo,zo,'*','LineWidth',10)
% plot3(xo,yo,-Q1zv,'*','LineWidth',10)
%-------------------------------------------
% yoo= linspace(-y_Room,y_Room,7);
% for j =1:length(yoo)
%     po=[xo,yo,zo];
%     [Xc,Yc,Zc] = cylinder(Ro);
%     Xc=Xc+po(1);
%     Yc=Yc+po(2);
%     Zc=1.7*Zc-z_Room;
%     surf(Xc,Yc,Zc)  
%     colormap summer
%     shading interp
%     hold on
% % end
% %Reflaction from obj
% plot3( [pt(1) pp(1)],[pt(2) pp(2)],[pt(3) pp(3)],'b','LineWidth',1.5)
% plot3( [pp(1) pr(1)],[pp(2) pr(2)],[pp(3) pr(3)],'b','LineWidth',1.5)
% plot3( pp(1),pp(2),pp(3),'X')
% figure(1)
% text(0.4,0.4,['Pr_total_dBm=',num2str(Pr_total_dBm)],'Units','normalized','Fontsize',nFont)
% hold off
%Diffraction from obj
% plot3( [pt(1) xo],[yt yo],[zt zo],'--','LineWidth',0.1)
% plot3( [xo x_Room],[yo yr],[zo zr],'--','LineWidth',0.1)
% plot3( [xo vec(1)],[yo vec(2)],[zo vec(3)],'--','LineWidth',0.1)
% hold off
%=========================#end#========================================

