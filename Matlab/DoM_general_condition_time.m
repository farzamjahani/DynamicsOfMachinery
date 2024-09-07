clc; clear;
%% INITIAL DATA
step = 64; accuracy = 9; n = 0;
time_initial = 0; time_final = pi/5; time_interval = time_final/step; %[s]
time = time_initial : time_interval : time_final; %[s]
teta2_initial = 0; %[rad]
w2 = 10; %[rad/s]
omega2 = [0 0 w2]; omega4 = [0 0 0]; %[rad/s]
alpha2 = [0 0 0]; alpha4 = [0 0 0]; %[rad/s^2]
% teta2 = teta initial + time*w2
teta2 = teta2_initial + time*w2; %[rad]
O2A = 0.15; AB = 0.25; BC = 0.25; AC = 0.5; BD = 0.28; O6D = 0.28; O2O6 = 0.25; %[m]
xO2 = 0; yO2 = 0; xO6 = 0; yO6 = 0.25; yC = 0; %[m]
MVx_A = zeros(1,step+1); MVx_B = zeros(1,step+1); MVx_C = zeros(1,step+1); MVx_D = zeros(1,step+1); %[m/s]
MVy_A = zeros(1,step+1); MVy_B = zeros(1,step+1); MVy_C = zeros(1,step+1); MVy_D = zeros(1,step+1); %[m/s]
MAx_A = zeros(1,step+1); MAx_B = zeros(1,step+1); MAx_C = zeros(1,step+1); MAx_D = zeros(1,step+1); %[m/s^2]
MAy_A = zeros(1,step+1); MAy_B = zeros(1,step+1); MAy_C = zeros(1,step+1); MAy_D = zeros(1,step+1); %[m/s^2]
OMEGA_2 = zeros(1,step+1); OMEGA_3 = zeros(1,step+1); OMEGA_4 = zeros(1,step+1); OMEGA_5 = zeros(1,step+1); OMEGA_6 = zeros(1,step+1); %[rad/s]
ALPHA_2 = zeros(1,step+1); ALPHA_3 = zeros(1,step+1); ALPHA_4 = zeros(1,step+1); ALPHA_5 = zeros(1,step+1); ALPHA_6 = zeros(1,step+1); %[rad/s^2]
R_F12 = zeros(3,step+1); R_F14 = zeros(3,step+1); R_F16 = zeros(3,step+1); R_F23 = zeros(3,step+1); R_F34 = zeros(3,step+1); R_F35 = zeros(3,step+1); R_F56 = zeros(3,step+1); %[N]
R_M2 = zeros(3,step+1); %[N.m]
TETA2 = 0 : (2*pi)/step : 2*pi; %[rad]
% neglecting friction in system
gravity_constant = 9.8; %[N/kg]
g = [0 -gravity_constant 0]; %[N/kg]
% linear density > 1 [g/mm] = 1 [kg/m]
LD = 1; %[kg/m]
% links mass
mass2 = LD*O2A; mass3 = LD*AC; mass4 = 0.5; mass5 = LD*BD; mass6 = LD*O6D; %[kg]
% Moment of inertia in center of mass
I2 = 1/12*mass2*O2A^2; I3 = 1/12*mass3*AC^2; I5 = 1/12*mass5*BD^2; I6 = 1/12*mass2*O2A^2; % I4 = ?
%% LOOP
for time = time_initial : time_interval : time_final %[s]
    n = n+1;
    teta2 = teta2_initial + time*w2; %[rad]
    %% POSITION ANALYSIS (p,r)
pO2 = [xO2 yO2 0];
pO6 = [xO6 yO6 0];
xA = O2A*cos(teta2); yA = O2A*sin(teta2);
pA = [xA yA 0];
% position C
syms O2C;
ec_O2C = (yA^2 + (xA+O2C)^2 == AC^2);
O2C = solve(ec_O2C,O2C); O2C = vpa(O2C,accuracy);
% O2C > 0
if O2C(1) < 0
    O2C(1) = [];
else
    O2C(2) = [];
end
xC = -O2C;
pC = [xC yC 0];
% position B , AB = BC , B is in the midle
xB = (xA+xC)/2; yB = (yA+yC)/2;
pB = [xB yB 0];
% position D
syms xD yD;
ec_xDyD = [ (xB-xD)^2 + (yB-yD)^2 == BD^2 ; (xO6-xD)^2 + (yO6-yD)^2 == O6D^2 ];
sol_xDyD = solve(ec_xDyD,xD,yD);
xD = sol_xDyD.xD; xD = vpa(xD,accuracy);
yD = sol_xDyD.yD; yD = vpa(yD,accuracy);
% xD min , yD max
if xD(1) >= xD(2)
    xD(1) = [];
else
    xD(2) = [];
end
if yD(1) <= yD(2)
    yD(1) = [];
else
    yD(2) = [];
end
pD = [xD yD 0];
rO2A = pA - pO2; rAC = pC - pA; rAB = pB - pA; rBD = pD - pB; rO6D = pD - pO6;
pO2 = vpa(pO2,accuracy); pO6 = vpa(pO6,accuracy); pA = vpa(pA,accuracy); pB = vpa(pB,accuracy); pC = vpa(pC,accuracy); pD = vpa(pD,accuracy);
rO2A = vpa(rO2A,accuracy); rAC = vpa(rAC,accuracy); rAB = vpa(rAB,accuracy); rBD = vpa(rBD,accuracy); rO6D = vpa(rO6D,accuracy);
    %% VELOCITY ANALYSIS (v)
% vA = cross(omega2,rO2A)
vA = cross(omega2,rO2A);
% vC = vA + vC/A , vC/A = cross(omega3,rAC)
syms vCx w3
vC = [vCx 0 0];
omega3 = [0 0 w3];
ec_vCw3 = (vC == vA + cross(omega3,rAC));
sol_vCw3 = solve(ec_vCw3,vCx,w3);
vC = [sol_vCw3.vCx 0 0];
omega3 = [0 0 sol_vCw3.w3];
% vB = vA + vB/A , vB/A = cross(omega3,rAB)
vB = vA + cross(omega3,rAB);
% vD = vB + vD/B , vD = cross(omega6,rO6D) , vD/B = cross(omega5,rBD)
syms w5 w6
omega5 = [0 0 w5];
omega6 = [0 0 w6];
ec_w5w6 = (cross(omega6,rO6D) == vB + cross(omega5,rBD));
sol_w5w6 = solve(ec_w5w6,w5,w6);
omega5 = [0 0 sol_w5w6.w5];
omega6 = [0 0 sol_w5w6.w6];
vD = cross(omega6,rO6D);
% Magnitude of Vector
mvA = (vA(1)^2 + vA(2)^2 + vA(3)^2)^(1/2);
mvB = (vB(1)^2 + vB(2)^2 + vB(3)^2)^(1/2);
mvC = (vC(1)^2 + vC(2)^2 + vC(3)^2)^(1/2);
mvD = (vD(1)^2 + vD(2)^2 + vD(3)^2)^(1/2);
vA = vpa(vA,accuracy); vB = vpa(vB,accuracy); vC = vpa(vC,accuracy); vD = vpa(vD,accuracy);
%%%%
    %% ACCELERATION ANALYSIS (a) aX = aXn + aXt
% aA = cross(omega2,cross(omega2,rO2A)) + cross(alpha2,rO2A)
aA = cross(omega2,cross(omega2,rO2A)) + cross(alpha2,rO2A);
% aC = aA + aA/C , aCn = 0 , aA/C = cross(omega3,cross(omega3,rAC)) + cross(alpha3,rAC)
syms aCx a3
aC = [aCx 0 0];
alpha3 = [0 0 a3];
ec_aCxa3 = (aC == aA + cross(omega3,cross(omega3,rAC)) + cross(alpha3,rAC));
sol_aCxa3 = solve(ec_aCxa3,aCx,a3);
aC = [sol_aCxa3.aCx 0 0];
alpha3 = [0 0 sol_aCxa3.a3];
% aB = aA + aB/A , aB/A = cross(omega3,cross(omega3,rAB)) + cross(alpha3,rAB)
aB = aA + cross(omega3,cross(omega3,rAB)) + cross(alpha3,rAB);
% aD = aB + aD/B , aD = cross(omega6,cross(omega6,rO6D)) + cross(alpha6,rO6D) , aD/B = cross(omega5,cross(omega5,rBD)) + cross(alpha5,rBD)
syms a5 a6
alpha5 = [0 0 a5];
alpha6 = [0 0 a6];
ec_a5a6 = (aB + cross(omega5,cross(omega5,rBD)) + cross(alpha5,rBD) == cross(omega6,cross(omega6,rO6D)) + cross(alpha6,rO6D));
sol_a5a6 = solve(ec_a5a6,a5,a6);
alpha5 = [0 0 sol_a5a6.a5];
alpha6 = [0 0 sol_a5a6.a6];
aD = cross(omega6,cross(omega6,rO6D)) + cross(alpha6,rO6D);
% Magnitude of Vector
maA = (aA(1)^2 + aA(2)^2 + aA(3)^2)^(1/2);
maB = (aB(1)^2 + aB(2)^2 + aB(3)^2)^(1/2);
maC = (aC(1)^2 + aC(2)^2 + aC(3)^2)^(1/2);
maD = (aD(1)^2 + aD(2)^2 + aD(3)^2)^(1/2);
aA = vpa(aA,accuracy); aB = vpa(aB,accuracy); aC = vpa(aC,accuracy); aD = vpa(aD,accuracy);
%%%%%%%%%%
    %% FORCE ANALYSIS (F,M)
rO2Amid = rO2A/2; rACmid = rAC/2; rBDmid = rBD/2; rO6Dmid = rO6D/2;
rc2 = rO2A/2; rc3 = pB; rc4 = pC; rc5 = (pB+pD)/2; rc6 = rO6D/2;
% acceleration in center of mass
ac2 = cross(omega2,cross(omega2,rO2Amid)) + cross(alpha2,rO2Amid);
ac3 = aA + cross(omega3,cross(omega3,rACmid)) + cross(alpha3,rACmid);
ac4 = aC;
ac5 = aB + cross(omega5,cross(omega5,rBDmid)) + cross(alpha5,rBDmid);
ac6 = cross(omega6,cross(omega6,rO6Dmid)) + cross(alpha6,rO6Dmid);
ac2 = vpa(ac2,12); ac3 = vpa(ac3,12); ac4 = vpa(ac4,12); ac5 = vpa(ac5,12); ac6 = vpa(ac6,12);
% all links
syms I4 F12x F12y M2z F14y F16x F16y F23x F23y F34x F34y F35x F35y F56x F56y
F12 = [F12x F12y 0]; F14 = [0 F14y 0]; F16 = [F16x F16y 0]; F23 = [F23x F23y 0]; F34 = [F34x F34y 0]; F35 = [F35x F35y 0]; F56 = [F56x F56y 0]; %[N]
M2 = [0 0 M2z]; %[N.m]
% Fab = -Fba
mg2 = mass2*g; mg3 = mass3*g; mg4 = mass4*g; mg5 = mass5*g; mg6 = mass6*g; %[N]
mac2 = mass2*ac2; mac3 = mass3*ac3; mac4 = mass4*ac4; mac5 = mass5*ac5; mac6 = mass6*ac6; %[N]
Ialpha2 = I2*alpha2; Ialpha3 = I3*alpha3; Ialpha4 = I4*alpha4; Ialpha5 = I5*alpha5; Ialpha6 = I6*alpha6; %[N.m]
% LINK (2)
ec_link2_F = (F12 + (-F23) + mg2 == mac2);
ec_link2_M = (cross(-rO2Amid,F12) + cross(rO2Amid,(-F23)) + cross([0 0 0],mg2) + M2 == Ialpha2);
% LINK (3)
ec_link3_F = (F23 + (-F34) + (-F35) + mg3 == mac3);
ec_link3_M = (cross(-rACmid,F23) + cross(rACmid,-F34) + cross([0 0 0],-F35) + cross([0 0 0],mg3) == Ialpha3);
% LINK (4)
ec_link4_F = (F34 + F14 + mg4 == mac4);
ec_link4_M = (cross([0 0 0],F34) + cross([0 0 0],F14) + cross([0 0 0],mg4) == Ialpha4);
% LINK (5)
ec_link5_F = (F35 + (-F56) + mg5 == mac5);
ec_link5_M = (cross(-rBDmid,F35) + cross(rBDmid,-F56) + cross([0 0 0],mg5) == Ialpha5);
% LINK (6)
ec_link6_F = (F16 + F56 + mg6 == mac6);
ec_link6_M = (cross(-rO6Dmid,F16) + cross(rO6Dmid,F56) + cross([0 0 0],mg6) == Ialpha6);
    %% FORCE RESULT
force_result = solve([ec_link2_F,ec_link2_M,ec_link3_F,ec_link3_M,ec_link4_F,ec_link4_M,ec_link5_F,ec_link5_M,ec_link6_F,ec_link6_M],[F12x F12y M2z F14y F16x F16y F23x F23y F34x F34y F35x F35y F56x F56y]);
F12x = force_result.F12x; F12y = force_result.F12y; M2z = force_result.M2z; F14y = force_result.F14y; F16x = force_result.F16x; F16y = force_result.F16y; F23x = force_result.F23x; F23y = force_result.F23y; F34x = force_result.F34x; F34y = force_result.F34y; F35x = force_result.F35x; F35y = force_result.F35y; F56x = force_result.F56x; F56y = force_result.F56y;
F12 = [F12x F12y 0]; F14 = [0 F14y 0]; F16 = [F16x F16y 0]; F23 = [F23x F23y 0]; F34 = [F34x F34y 0]; F35 = [F35x F35y 0]; F56 = [F56x F56y 0]; %[N]
M2 = [0 0 M2z]; %[N.m]
    %% PLOT ANIMATION
P = subplot(2,2,1);
V = subplot(2,2,2);
A = subplot(2,2,4);
% POSITION
POSITION = plot(P,[xO2 xA],[yO2 yA],[xA xC],[yA yC],[xB xD],[yB yD],[xO6 xD],[yO6 yD]);
axis(P,[-0.7 0.2 -0.2 0.4]);
grid(P,'on'); grid(P,'minor');
xlabel(P,'X (m)'); ylabel(P,'Y (m)');
title(P,'POSITION ANALYSIS');
text(P,xO2,yO2,'O2'); text(P,xO6,yO6,'O6'); text(P,xC,yC,'C'); text(P,xA,yA,'A'); text(P,xB,yB,'B'); text(P,xD,yD,'D');
POSITION(1).LineWidth = 1.6; POSITION(2).LineWidth = 1.6; POSITION(3).LineWidth = 1.6; POSITION(4).LineWidth = 1.6;
% VELOCITY
VELOCITY = plot(V,[0 vA(1)],[0 vA(2)],[0 vB(1)],[0 vB(2)],[0 vC(1)],[0 vC(2)],[0 vD(1)],[0 vD(2)]);
axis(V,[-1.75 1.75 -1.75 1.75]);
grid(V,'on'); grid(V,'minor');
xlabel(V,'Vx (m/s)'); ylabel(V,'Vy (m/s)');
title(V,'VELOCITY DIAGRAM');
text(V,0,0,'Ov'); text(V,vA(1),vA(2),'A'); text(V,vB(1),vB(2),'B'); text(V,vC(1),vC(2),'C'); text(V,vD(1),vD(2),'D');
VELOCITY(1).LineWidth = 1.6; VELOCITY(2).LineWidth = 1.6; VELOCITY(3).LineWidth = 1.6; VELOCITY(4).LineWidth = 1.6;
% ACCELERATION
ACCELERATION = plot(A,[0 aA(1)],[0 aA(2)],[0 aB(1)],[0 aB(2)],[0 aC(1)],[0 aC(2)],[0 aD(1)],[0 aD(2)]);
axis(A,[-17 22 -20 17]);
grid(A,'on'); grid(A,'minor');
xlabel(A,'Ax (m/s^2)'); ylabel(A,'Ay (m/s^2)');
title(A,'ACCELERATION DIAGRAM');
text(A,0,0,'Oa'); text(A,aA(1),aA(2),'A'); text(A,aB(1),aB(2),'B'); text(A,aC(1),aC(2),'C'); text(A,aD(1),aD(2),'D');
ACCELERATION(1).LineWidth = 1.6; ACCELERATION(2).LineWidth = 1.6; ACCELERATION(3).LineWidth = 1.6; ACCELERATION(4).LineWidth = 1.6;
    %% PAUSE
pause(0.0001) %(time_interval)
    %% SAVE VELOCITY & ACCELERATION & OMEGA & ALPHA
MVx_A(n) = vA(1); MVx_B(n) = vB(1); MVx_C(n) = vC(1); MVx_D(n) = vD(1); %[m/s]
MVy_A(n) = vA(2); MVy_B(n) = vB(2); MVy_C(n) = vC(2); MVy_D(n) = vD(2); %[m/s]
MAx_A(n) = aA(1); MAx_B(n) = aB(1); MAx_C(n) = aC(1); MAx_D(n) = aD(1); %[m/s^2]
MAy_A(n) = aA(2); MAy_B(n) = aB(2); MAy_C(n) = aC(2); MAy_D(n) = aD(2); %[m/s^2]
OMEGA_2(n) = omega2(3); OMEGA_3(n) = omega3(3); OMEGA_4(n) = omega4(3); OMEGA_5(n) = omega5(3); OMEGA_6(n) = omega6(3); %[rad/s]
ALPHA_2(n) = alpha2(3); ALPHA_3(n) = alpha3(3); ALPHA_4(n) = alpha4(3); ALPHA_5(n) = alpha5(3); ALPHA_6(n) = alpha6(3); %[rad/s^2]
R_F12(:,n) = F12; R_F14(:,n) = F14; R_F16(:,n) = F16; R_F23(:,n) = F23; R_F34(:,n) = F34; R_F35(:,n) = F35; R_F56(:,n) = F56; %[N]
R_M2(:,n) = M2; %[N.m]
part = n
end
time = time_initial : time_interval : time_final;
%% PLOT RESULT (V,A)
figure
Vrx = subplot(2,3,1);
Vry = subplot(2,3,2);
OMEGArz = subplot(2,3,3);
Arx = subplot(2,3,4);
Ary = subplot(2,3,5);
ALPHArz = subplot(2,3,6);
% VELOCITY X
VELOCITY_Rx = plot(Vrx,time,MVx_A,time,MVx_B,time,MVx_C,time,MVx_D);
grid(Vrx,'on'); grid(Vrx,'minor');
xlabel(Vrx,'time (s)'); ylabel(Vrx,'Vx (m/s)');
VX=legend(Vrx,'A','B','C','D'); title(VX,'point');
title(Vrx,'VELOCITY IN X DIRECTION');
VELOCITY_Rx(1).LineWidth = 1; VELOCITY_Rx(2).LineWidth = 1; VELOCITY_Rx(3).LineWidth = 1; VELOCITY_Rx(4).LineWidth = 1;
% VELOCITY Y
VELOCITY_Ry = plot(Vry,time,MVy_A,time,MVy_B,time,MVy_C,time,MVy_D);
grid(Vry,'on'); grid(Vry,'minor');
xlabel(Vry,'time (s)'); ylabel(Vry,'Vy (m/s)');
VY=legend(Vry,'A','B','C','D'); title(VY,'point');
title(Vry,'VELOCITY IN Y DIRECTION');
VELOCITY_Ry(1).LineWidth = 1; VELOCITY_Ry(2).LineWidth = 1; VELOCITY_Ry(3).LineWidth = 1; VELOCITY_Ry(4).LineWidth = 1;
% OMEGA Z
OMEGA_Rz = plot(OMEGArz,time,OMEGA_2,time,OMEGA_3,time,OMEGA_4,time,OMEGA_5,time,OMEGA_6);
grid(OMEGArz,'on'); grid(OMEGArz,'minor');
xlabel(OMEGArz,'time (s)'); ylabel(OMEGArz,'OMEGAz (rad/s)');
OZ=legend(OMEGArz,'2','3','4','5','6'); title(OZ,'link');
title(OMEGArz,'OMEGA IN Z DIRECTION');
OMEGA_Rz(1).LineWidth = 1; OMEGA_Rz(2).LineWidth = 1; OMEGA_Rz(3).LineWidth = 1; OMEGA_Rz(4).LineWidth = 1; OMEGA_Rz(5).LineWidth = 1;
% ACCELERATION X
ACCELERATION_Rx = plot(Arx,time,MAx_A,time,MAx_B,time,MAx_C,time,MAx_D);
grid(Arx,'on'); grid(Arx,'minor');
xlabel(Arx,'time (s)'); ylabel(Arx,'Ax (m/s^2)');
AX=legend(Arx,'A','B','C','D'); title(AX,'point');
title(Arx,'ACCELERATION IN X DIRECTION');
ACCELERATION_Rx(1).LineWidth = 1; ACCELERATION_Rx(2).LineWidth = 1; ACCELERATION_Rx(3).LineWidth = 1; ACCELERATION_Rx(4).LineWidth = 1;
% ACCELERATION Y
ACCELERATION_Ry = plot(Ary,time,MAy_A,time,MAy_B,time,MAy_C,time,MAy_D);
grid(Ary,'on'); grid(Ary,'minor');
xlabel(Ary,'time (s)'); ylabel(Ary,'Ay (m/s^2)');
AY=legend(Ary,'A','B','C','D'); title(AY,'point');
title(Ary,'ACCELERATION IN Y DIRECTION');
ACCELERATION_Ry(1).LineWidth = 1; ACCELERATION_Ry(2).LineWidth = 1; ACCELERATION_Ry(3).LineWidth = 1; ACCELERATION_Ry(4).LineWidth = 1;
% ALPHA Z
ALPHA_Rz = plot(ALPHArz,time,ALPHA_2,time,ALPHA_3,time,ALPHA_4,time,ALPHA_5,time,ALPHA_6);
grid(ALPHArz,'on'); grid(ALPHArz,'minor');
xlabel(ALPHArz,'T (rad)'); ylabel(ALPHArz,'ALPHAz (rad/s^2)');
AZ=legend(ALPHArz,'2','3','4','5','6'); title(AZ,'link');
title(ALPHArz,'ALPHA IN Z DIRECTION');
ALPHA_Rz(1).LineWidth = 1; ALPHA_Rz(2).LineWidth = 1; ALPHA_Rz(3).LineWidth = 1; ALPHA_Rz(4).LineWidth = 1; ALPHA_Rz(5).LineWidth = 1;
%% PLOT RESULT (F) (x,y)
figure
Fx = subplot(2,1,1);
Fy = subplot(2,1,2);
% F x
F_Rx = plot(Fx,time,R_F12(1,:),time,R_F14(1,:),time,R_F16(1,:),time,R_F23(1,:),time,R_F34(1,:),time,R_F35(1,:),time,R_F56(1,:));
grid(Fx,'on'); grid(Fx,'minor');
xlabel(Fx,'time (s)'); ylabel(Fx,'force x (N)');
FX=legend(Fx,'F12','F14','F16','F23','F34','F35','F56'); title(FX,'force');
title(Fx,'FORCE IN X DIRECTION');
F_Rx(1).LineWidth = 1; F_Rx(2).LineWidth = 1; F_Rx(3).LineWidth = 1; F_Rx(4).LineWidth = 1; F_Rx(5).LineWidth = 1; F_Rx(6).LineWidth = 1; F_Rx(7).LineWidth = 1;
% F y
F_Ry = plot(Fy,time,R_F12(2,:),time,R_F14(2,:),time,R_F16(2,:),time,R_F23(2,:),time,R_F34(2,:),time,R_F35(2,:),time,R_F56(2,:));
grid(Fy,'on'); grid(Fy,'minor');
xlabel(Fy,'time (s)'); ylabel(Fy,'force y (N)');
FY=legend(Fy,'F12','F14','F16','F23','F34','F35','F56'); title(FY,'force');
title(Fy,'FORCE IN Y DIRECTION');
F_Ry(1).LineWidth = 1; F_Ry(2).LineWidth = 1; F_Ry(3).LineWidth = 1; F_Ry(4).LineWidth = 1; F_Ry(5).LineWidth = 1; F_Ry(6).LineWidth = 1; F_Ry(7).LineWidth = 1;
%% PLOT RESULT (M) (z)
figure
Mz = subplot(1,1,1);
M_Ry = plot(Mz,time,R_M2(3,:));
grid(Mz,'on'); grid(Mz,'minor');
xlabel(Mz,'time (s)'); ylabel(Mz,'moment z (N.m)');
MZ=legend(Mz,'M2'); title(MZ,'moment');
title(Fy,'MOMENT IN Z DIRECTION');
M_Ry(1).LineWidth = 1;