% Project 2 - EAE 135
% Team 20
% DUE: 02/24/26



% Given Parameters
F = 726000;                     % max thrust of first stage [N]
rocket_diam = 1.28;             % rocket diameter [m]
prop_diam = 1.18;               % propellant diameter [m]
AoA = 20;                       % angle of attack [deg]
length = 5*rocket_diam;         % rocket length [m]
thickness_oneLayer = 0.050 / 8;        % thickness of one layer of carbon epoxy (m)
cpcg = rocket_diam;             % distance btwn cp-cg [m]
P1_net = 600000;                % net axial force on O [N]
Weight = 8500;                  % weight of rocket [kg]
lift = 1400000;                 % lift force [N]
FoS = 1.25;                     % spacecraft safety factor

% AS4/Epoxy properties
rho = 1522.3948;                            % density [kg/m^3]
E11 = 148e3;                                % Axial Young's Modulus [Mpa]
E22 = 10.5e3;                               % Transverse Young's Modulus [Mpa]
v12 = 0.30;                                 % Major Poisson's ratio
v21 = v12*E22/E11;                          % Minor Poisson's ratio
G12 = 5.61e3;                               % Shear modulus [MPa]

% Strengths for failure analysis
StrengthTensile_0deg = 2137;                % MPa - Axial Tensile Strength (0 deg ply)
StrengthCompressive_0deg = -184 * 6.89476;  % -184 ksi to MPa (0 deg ply)
StrengthTensile_90deg = 53.4;               % Mpa
StrengthCompressive_90deg = -24.4 * 6.89476;% -24.4 ksi to MPa (90 deg ply)
StrengthTensile_45deg = 53.4;               % Mpa
StrengthCompressive_45deg = -24.4 * 6.89476;% -24.4 ksi to MPa (+/-45 deg ply)

% Propellant Properties (HTPB, Aged 185 Days)



% Layup [0/+/-45/90]s (8 layers of CFRP)
% ----------------- [0 deg (r=1.28/2 = 0.64m)]
% ----------------- [+45 deg]
% ----------------- [-45 deg]
% ----------------- [90 deg]
% ----------------- [90 deg]
% ----------------- [-45 deg]
% ----------------- [+45 deg]
% ----------------- [0 deg]
% \\\\\\\\\\\\\\\\\ [CORE (r=1.18/2 = 0.59m)]
% \\\\\\  ^  \\\\\\
% \\\\\\  |  \\\\\\
% +++++++++++++++++ [Centerline (r=0)]
% Q Matrix (Same for all layers)
Q = [ E11/(1-v12*v21), v21*E11/(1-v12*v21), 0;
      v12*E22/(1-v12*v21), E22/(1-v12*v21), 0;
      0,                    0,            G12]

% Find Qbar for each layer
beta_0deg = deg2rad(0);
beta_45deg = deg2rad(45);
beta_90deg = deg2rad(90);

function T = calcTMatrixFromBeta(beta)
    T = [cos(beta)^2, sin(beta)^2, 2*sin(beta)*cos(beta);
         sin(beta)^2, cos(beta)^2, -2*sin(beta)*cos(beta);
         -1*sin(beta)*cos(beta), sin(beta)*cos(beta), cos(beta)^2 - sin(beta)^2];
end

T_0deg = calcTMatrixFromBeta(beta_0deg);
T_pos45deg = calcTMatrixFromBeta(beta_45deg);
T_neg45deg = calcTMatrixFromBeta(-beta_45deg);
T_90deg = calcTMatrixFromBeta(beta_90deg);

R = [1 0 0; 0 1 0; 0 0 2];

function Qbar = calcQbar(Q, T, R)
    T_inv = inv(T);
    R_inv = inv(R);
    Qbar = T_inv * Q * R * T * R_inv;
end

Qbar_0deg = calcQbar(Q, T_0deg, R)
Qbar_pos45deg = calcQbar(Q, T_pos45deg, R)
Qbar_neg45deg = calcQbar(Q, T_neg45deg, R)
Qbar_90deg = calcQbar(Q, T_90deg, R)

% QBAR VALUE CHECK
% Q0deg (1,1) = Q90deg (2,2)   - PASS
% Q45deg (1,1) = Q45deg (2,2)  - PASS
% Q45deg (3,3) > Q0deg (3,3)   - PASS


Sbar_0deg = inv(Qbar_0deg);
Sbar_pos45deg = inv(Qbar_pos45deg);
Sbar_neg45deg = inv(Qbar_neg45deg);
Sbar_90deg = inv(Qbar_90deg);

% Extract projected Young's modulus from Sbar (Ex1x1)


% Iterate to find H33_C and Stress per layer