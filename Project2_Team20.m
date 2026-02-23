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
rho = 1522.3948;            % density (kg/m^3)
E11 = 21.5*6.895e9;         % longitudinal modulus (Pa)
E22 = 1.46*6.895e9;         % transverse modulus (Pa)
v12 = 0.30;                 % major Poisson's ratio
G12 = 0.81*6.895e9;         % shear modulus (Pa)
tens_ax = 310*6.895e6;      % axial tensile strength (Pa)
comp_ax = -184*6.895e6;     % axial compressive strength (Pa)
tens_trav = 7.75*6.895e6;   % transverse tensile strength (Pa)
comp_trav = -24.4*6.895e6;  % transverse compressive strength (Pa)
strain_tens_ax = 0.014;     % axial tensile failure strain
strain_comp_ax = -0.01;     % axial compressive failure strain

% Axial Stress/Displacement
A = pi*(rk_d/2)^2 - pi*(prop_d/2)^2;  % cross-sectional area of the case (m^2)

E0  = E11;
E90 = E22;
E45 = 1/(1/E11 + 1/E22);
Eeff = (2*E0 + 4*E45 + 2*E90)/8;        % Simple effective modulus for [0/+45/-45/90]s *idk maybe*

P_tot = P - lift*sin(aoa);
u1_axial = @(x1) (P_tot/(Eeff*A))*(x1 - l);

u_left = u1_axial(0);       % max magnitude
u_right = u1_axial(l);      % fixed end

x1 = linspace(0, l, 200);
u_axial = u1_axial(x1);

figure;
plot(x1, u_axial, 'LineWidth', 1.5); grid on;
xlabel('x_1 (m)'); ylabel('u_{1,axial} (m)');
title('Axial displacement along beam');