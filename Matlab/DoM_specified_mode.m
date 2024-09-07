clc; clear;
%% SYSTEM INITIAL DATA
O2A = 0.15; AB = 0.25; BC = 0.25; AC = 0.5; BD = 0.28; O6D = 0.28; O2O6 = 0.25; %[m]
xO2 = 0; yO2 = 0; xO6 = 0; yO6 = 0.25; yC = 0; %[m]
teta2 = pi/4; %[rad]
w2 = 10; %[rad/s]
omega2 = [0 0 w2]; omega4 = [0 0 0]; %[rad/s]
alpha2 = [0 0 0]; alpha4 = [0 0 0]; %[rad/s^2]
%% POSITION ANALYSIS (p,r)
pO2 = [xO2 yO2 0];
pO6 = [xO6 yO6 0];
xA = O2A*cos(teta2); yA = O2A*sin(teta2);
pA = [xA yA 0];
% position C
syms O2C;
ec_O2C = (yA^2 + (xA+O2C)^2 == AC^2);
O2C = solve(ec_O2C,O2C);
O2C = vpa(O2C,12);
% O2c > 0
O2C(1) = [];
xC = -O2C;
pC = [xC yC 0];
% position B , AB = BC , B is in the midle
xB = (xA+xC)/2; yB = (yA+yC)/2;
pB = [xB yB 0];
% position D
syms xD yD;
ec_xDyD = [ (xB-xD)^2 + (yB-yD)^2 == BD^2 ; (xO6-xD)^2 + (yO6-yD)^2 == O6D^2 ];
sol_xDyD = solve(ec_xDyD,xD,yD);
xD = sol_xDyD.xD;
yD = sol_xDyD.yD;
% xD < 0 , yD > yB
xD(2) = [];
yD(2) = [];
pD = [xD yD 0];
rO2A = pA - pO2; rAC = pC - pA; rAB = pB - pA; rBD = pD - pB; rO6D = pD - pO6;
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
%% FORCE ANALYSIS (F,M)
% neglecting friction in system
gravity_constant = 9.8; %[N/kg]
g = [0 -gravity_constant 0]; %[N/kg]
% linear density > 1 [g/mm] = 1 [kg/m]
LD = 1; %[kg/m]
% links mass
mass2 = LD*O2A; mass3 = LD*AC; mass4 = 0.5; mass5 = LD*BD; mass6 = LD*O6D; %[kg]
% moment of inertia in center of mass
I2 = 1/12*mass2*O2A^2; I3 = 1/12*mass3*AC^2; I5 = 1/12*mass5*BD^2; I6 = 1/12*mass2*O2A^2; % I4 = ?
rO2Amid = rO2A/2; rACmid = rAC/2; rBDmid = rBD/2; rO6Dmid = rO6D/2;
rc2 = rO2A/2; rc3 = pB; rc4 = pC; rc5 = (pB+pD)/2; rc6 = rO6D/2;
% acceleration in center of mass
ac2 = cross(omega2,cross(omega2,rO2Amid)) + cross(alpha2,rO2Amid);
ac3 = aA + cross(omega3,cross(omega3,rACmid)) + cross(alpha3,rACmid);
ac4 = aC;
ac5 = aB + cross(omega5,cross(omega5,rBDmid)) + cross(alpha5,rBDmid);
ac6 = cross(omega6,cross(omega6,rO6Dmid)) + cross(alpha6,rO6Dmid);
%% ALL LINKS
syms I4 F12x F12y M2z F14y F16x F16y F23x F23y F34x F34y F35x F35y F56x F56y
F12 = [F12x F12y 0]; F14 = [0 F14y 0]; F16 = [F16x F16y 0]; F23 = [F23x F23y 0]; F34 = [F34x F34y 0]; F35 = [F35x F35y 0]; F56 = [F56x F56y 0]; %[N]
M2 = [0 0 M2z]; %[N.m]
% Fab = -Fba
mg2 = mass2*g; mg3 = mass3*g; mg4 = mass4*g; mg5 = mass5*g; mg6 = mass6*g; %[N]
mac2 = mass2*ac2; mac3 = mass3*ac3; mac4 = mass4*ac4; mac5 = mass5*ac5; mac6 = mass6*ac6; %[N]
Ialpha2 = I2*alpha2; Ialpha3 = I3*alpha3; Ialpha4 = I4*alpha4; Ialpha5 = I5*alpha5; Ialpha6 = I6*alpha6; %[N.m]
%% LINK (2)
ec_link2_F = (F12 + (-F23) + mg2 == mac2);
ec_link2_M = (cross(-rO2Amid,F12) + cross(rO2Amid,(-F23)) + cross([0 0 0],mg2) + M2 == Ialpha2);
%% LINK (3)
ec_link3_F = (F23 + (-F34) + (-F35) + mg3 == mac3);
ec_link3_M = (cross(-rACmid,F23) + cross(rACmid,-F34) + cross([0 0 0],-F35) + cross([0 0 0],mg3) == Ialpha3);
%% LINK (4)
ec_link4_F = (F34 + F14 + mg4 == mac4);
ec_link4_M = (cross([0 0 0],F34) + cross([0 0 0],F14) + cross([0 0 0],mg4) == Ialpha4);
%% LINK (5)
ec_link5_F = (F35 + (-F56) + mg5 == mac5);
ec_link5_M = (cross(-rBDmid,F35) + cross(rBDmid,-F56) + cross([0 0 0],mg5) == Ialpha5);
%% LINK (6)
ec_link6_F = (F16 + F56 + mg6 == mac6);
ec_link6_M = (cross(-rO6Dmid,F16) + cross(rO6Dmid,F56) + cross([0 0 0],mg6) == Ialpha6);
%% FORCE RESULT
force_result = solve([ec_link2_F,ec_link2_M,ec_link3_F,ec_link3_M,ec_link4_F,ec_link4_M,ec_link5_F,ec_link5_M,ec_link6_F,ec_link6_M],[F12x F12y M2z F14y F16x F16y F23x F23y F34x F34y F35x F35y F56x F56y]);
F12x = force_result.F12x; F12y = force_result.F12y; M2z = force_result.M2z; F14y = force_result.F14y; F16x = force_result.F16x; F16y = force_result.F16y; F23x = force_result.F23x; F23y = force_result.F23y; F34x = force_result.F34x; F34y = force_result.F34y; F35x = force_result.F35x; F35y = force_result.F35y; F56x = force_result.F56x; F56y = force_result.F56y;
F12 = [F12x F12y 0]; F14 = [0 F14y 0]; F16 = [F16x F16y 0]; F23 = [F23x F23y 0]; F34 = [F34x F34y 0]; F35 = [F35x F35y 0]; F56 = [F56x F56y 0]; %[N]
M2 = [0 0 M2z]; %[N.m]
%% RESULT PLOT
R_O2C = vpa(O2C,12) %[O2C]
R_xD = vpa(xD,12) , R_yD = vpa(xD,12) %[D]
R_pO2 = vpa(pO2,12) , R_pO6 = vpa(pO6,12) , R_pA = vpa(pA,12) , R_pB = vpa(pB,12) , R_pC = vpa(pC,12) , R_pD = vpa(pD,12) %[p]
R_rO2A = vpa(rO2A,12) , R_rAC = vpa(rAC,12) , R_rAB = vpa(rAB,12) , R_rBD = vpa(rBD,12) , R_rO6D = vpa(rO6D,12) %[r]
R_vA = vpa(vA,12) , R_vB = vpa(vB,12) , R_vC = vpa(vC,12) , R_vD = vpa(vD,12) %[v]
R_mvA = vpa(mvA,12) , R_mvB = vpa(mvB,12) , R_mvC = vpa(mvC,12) ,R_mvD = vpa(mvD,12) %[mv]
R_aA = vpa(aA,12) , R_aB = vpa(aB,12) , R_aC = vpa(aC,12) , R_aD = vpa(aD,12) %[a]
R_maA = vpa(maA,12) , R_maB = vpa(maB,12) , R_maC = vpa(maC,12) , R_maD = vpa(maD,12) %[ma]
R_ac2 = vpa(ac2,12) , R_ac3 = vpa(ac3,12) , R_ac4 = vpa(ac4,12) , R_ac5 = vpa(ac5,12) , R_ac6 = vpa(ac6,12) %[ac]
R_F12 = vpa(F12,12) , R_F14 = vpa(F14,12) , R_F16 = vpa(F16,12) , R_F23 = vpa(F23,12) , R_F34 = vpa(F34,12) , R_F35 = vpa(F35,12) , R_F56 = vpa(F56,12) %[F]
R_M2 = vpa(M2,12) %[M]


