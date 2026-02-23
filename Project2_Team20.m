% Project 2
% Team 20
% 2024-06-01

F = 726000;         % max thrust of first stage (N)
rk_d = 1.28;        % rocket diameter (m)
prop_d = 1.18;      % propellant diameter (m)
aoa = 20;           % angle of attack (deg)
l = 5*rk_d;         % rocket length (m)
t = 6.25e-3;        % thickness of one layer of carbon epoxy (m)
cpcg = rk_d;        % distance btwn cp-cg (m)
P = 600000;         % net axial force on O (N)
W = 8500;           % weight of rocket (kg)
lift = 1400000;     % lift force (N)
FoS = 1.25;         % safety factor

% AS4/Epoxy properties
rho = 1522.3948;        % density (kg/m^3)
E11 = 21.5*6.895e9;     % longitudinal modulus (Pa)
E22 = 1.46*6.895e9;     % transverse modulus (Pa)
v12 = 0.30;             % major Poisson's ratio
G12 = 0.81*6.895e9;     % shear modulus (Pa)
