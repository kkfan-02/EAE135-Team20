% Project 2 - EAE 135
% Team 20
% DUE: 02/24/26



% Given Parameters
F = 726000;                     % max thrust of first stage [N]
rocket_diam = 1.28;             % rocket diameter [m]
prop_diam = 1.18;               % propellant diameter [m]
AoA = 20;                       % angle of attack [deg]
beam_length = 5*rocket_diam;         % beam length [m]
thickness_oneLayer = 0.050 / 8;        % thickness of one layer of carbon epoxy (m)
cpcg = rocket_diam;             % distance btwn cp-cg [m]
P1_net = 600000;                % net axial force on O [N]
Weight = 8500;                  % weight of rocket [kg]
lift = 1400000;                 % lift force [N]
FoS = 1.25;                     % spacecraft safety factor
g = 9.81;                       % gravitational acceleration [m/s^2]

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
    0,                    0,            G12];

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
% Q0deg = Q                    - PASS
% Q0deg (1,1) = Q90deg (2,2)   - PASS
% Q45deg (1,1) = Q45deg (2,2)  - PASS
% Q45deg (3,3) > Q0deg (3,3)   - PASS

Sbar_0deg = inv(Qbar_0deg);
Sbar_pos45deg = inv(Qbar_pos45deg);
Sbar_neg45deg = inv(Qbar_neg45deg);
Sbar_90deg = inv(Qbar_90deg);

% Extract projected Young's modulus from Sbar (Ex1x1)
Ex1x1_0deg = 1 / (Sbar_0deg(1,1));          % Projected Young's Modulus for 0 deg layer
Ex1x1_pos45deg = 1 / (Sbar_pos45deg(1,1));  % Projected Young's Modulus for positive 45 deg layer
Ex1x1_neg45deg = 1 / (Sbar_neg45deg(1,1));  % Projected Young's Modulus for negative 45 deg layer
Ex1x1_90deg = 1 / (Sbar_90deg(1,1));        % Projected Young's Modulus for 90 deg layer

% ************************************************************************
% ************************** H33_C Loop **************************
n_layers = 8;
t = thickness_oneLayer;     % thickness of one ply [m]

% x2 boundaries for each layer, counting upward (inside to outside)
% x2_bounds(k)   = bottom (inner) face of layer k  -> x2_i
% x2_bounds(k+1) = top (outer) face of layer k     -> x2_(i+1)
% Layer 1 is the innermost ply (0 deg), layer 8 is the outermost ply (0 deg)
r_inner_wall = prop_diam / 2;                    % innermost radius (bottom of layer 1) [m]
x2_bounds = r_inner_wall + (0:n_layers) * t;     % [x2_1, x2_2, ..., x2_9]

% Layup ordered inside-out to match x2 counting upward:
% Layer 1 (innermost) -> Layer 8 (outermost)
ply_angles_deg = [0, 45, -45, 90, 90, -45, 45, 0];  % [0/+45/-45/90]s, inner to outer
Ex1x1_layers   = [Ex1x1_0deg, Ex1x1_pos45deg, Ex1x1_neg45deg, Ex1x1_90deg, ...
    Ex1x1_90deg, Ex1x1_neg45deg, Ex1x1_pos45deg, Ex1x1_0deg];

% --- Step 5: Iterate layer integrals to build H33_C ---
H33_C = 0;  % initialize bending stiffness [MPa·m^4]

for k = 1:n_layers
    x2_i   = x2_bounds(k);     % bottom (inner) face of layer k [m]
    x2_ip1 = x2_bounds(k+1);   % top (outer) face of layer k   [m]
    Ex     = Ex1x1_layers(k);  % projected Young's modulus of layer k [MPa]

    % Analytical result of: integral_0^2pi sin^2(theta) dtheta = pi
    % times: integral_{x2_i}^{x2_(i+1)} r^3 dr = (x2_(i+1)^4 - x2_i^4) / 4
    radial_integral = (x2_ip1^4 - x2_i^4) / 4;

    H33_C = H33_C + Ex * pi * radial_integral;
end

% H33_C is now the effective bending stiffness of the composite cylinder [MPa·m^4]
% (units: [MPa] * [m^4] = [N/m^2 * m^4] = [N·m^2])
disp(H33_C + " [MPa] * [m^4]")


% ************************************************************************
% ************** INSTRUCTIONS FOR NEXT ***********************************
% ************************************************************************
% M3 code probably not correct
% Refer to use_for_project2.pdf on canvas to find how to find u2 and u1 from here
% Also finish the accompanying video


% ***************** AXIAL PROBLEM *****************
% u1_total = u1_axial + u1_bending

% Axial Stiffness
S_axial = 0;

for k = 1:n_layers
    r_i = x2_bounds(k);         % inner radius of layer k [m]
    r_ip1 = x2_bounds(k+1);     % outer radius of layer k [m]
    Ex_MPa = Ex1x1_layers(k);   % projected Young's modulus of layer k [MPa]
    Ex_Pa = Ex_MPa * 1e6;       % convert MPa to Pa

    A_k = pi * (r_ip1^2 - r_i^2);  % cross-sectional area of layer k [m^2]
    S_axial = S_axial + Ex_Pa * A_k;  % sum axial stiffness contributions [N]
end

x1_vec = linspace(0, beam_length, 500);  % x1 positions along the beam [m]

% Axial Strain
u1prime_axial = P1_net / S_axial;  % du1/dx1 (dimensionless)

% Axial displacement with clamp @ x = rocket_length
u1_axial = u1prime_axial * (x1_vec - beam_length);  % [m]

disp("axial strain = " + u1prime_axial);
disp("u1_axial at x=0 (fixed end) = " + u1_axial(1) * 1000 + " mm");

% ***************** BENDING PROBLEM *****************

L_transverse = lift*cosd(AoA);     % [N]
W_transverse = Weight * g;         % [N]

H = H33_C * 1e6;                   % convert MPa·m^4 to N·m^2

% Build bending moment M3(x) piecewise
M3 = zeros(size(x1_vec));  % [N·m]
for i = 1:length(x1_vec)
    x = x1_vec(i);
    if x < cpcg
        M3(i) = L_transverse * x;
    else
        M3(i) = L_transverse * x - W_transverse * (x - cpcg);
    end
end

% Now solve for u2 using piecewise integrated forms with constants
% Unknown constants: C1, C2 (region1), C3, C4 (region 2)
% Equations
% (1) u2_2'(L)=0
% (2) u2_2(L)=0
% (3) u2_1'(cpcg)=u2_2'(cpcg)
% (4) u2_1(cpcg)=u2_2(cpcg)

% Region 1:
% u2_1'(x) = (L_transverse*x^2)/(2H) + C1
% u2_1(x)  = (L_transverse*x^3)/(6H) + C1*x + C2

% Region 2:
% u2_2'(x) = (L_transverse*x^2)/(2H) - (W_transverse*(x-cpcg)^2)/(2H) + C3
% u2_2(x)  = (L_transverse*x^3)/(6H) - (W_transverse*(x-cpcg)^3)/(6H) + C3*x + C4

A = zeros(4,4);
b = zeros(4,1);

% (1) u2_2'(L)=0:  (L*term) + C3 = 0  -> C3 coefficient 1
A(1,:) = [0 0 1 0];
b(1) = -(L_transverse*beam_length^2)/(2*H) + (W_transverse*(beam_length-cpcg)^2)/(2*H);

% (2) u2_2(L) = 0: C3*L + C4 = -(L_transverse*L^3) + (W_transverse*(L-cpcg)^3) / (6H)
A(2,:) = [0 0 beam_length 1];
b(2) = -(L_transverse*beam_length^3)/(6*H) + (W_transverse*(beam_length-cpcg)^3)/(6*H);

% (3) u2_1'(cpcg)=u2_2'(cpcg)
% Left:  (L*cpcg^2)/(2H) + C1
% Right: (L*cpcg^2)/(2H) - (W_transverse*(cpcg-cpcg)^2)/(2H) + C3 = (L*cpcg^2)/(2H) + C3 => C1 - C3 = 0
A(3,:) = [1 0 -1 0];
b(3) = 0;

% (4) u2_1(cpcg)=u2_2(cpcg)
% Left:  (L*cpcg^3)/(6H) + C1*cpcg + C2
% Right: (L*cpcg^3)/(6H) - (W_transverse*(cpcg-cpcg)^3)/(6H) + C3*cpcg + C4 = (L*cpcg^3)/(6H) + C3*cpcg + C4 => C1*cpcg + C2 - C3*cpcg - C4 = 0
A(4,:) = [cpcg 1 -cpcg -1];
b(4) = 0;

C = A\b;
C1 = C(1); C2 = C(2); C3 = C(3); C4 = C(4);

% Evaluate u2 and u2' at x1_vec
u2 = zeros(size(x1_vec));           % [m]
u2_prime = zeros(size(x1_vec));     % slope [dimensionless]

for i = 1: length(x1_vec)
    x = x1_vec(i);
    if x < cpcg
        u2_prime(i) = (L_transverse*x^2)/(2*H) + C1;
        u2(i) = (L_transverse*x^3)/(6*H) + C1*x + C2;
    else
        u2_prime(i) = (L_transverse*x^2)/(2*H) - (W_transverse*(x-cpcg)^2)/(2*H) + C3;
        u2(i) = (L_transverse*x^3)/(6*H) - (W_transverse*(x-cpcg)^3)/(6*H) + C3*x + C4;
    end
end

% BC checks
disp("u2_prime at x=0 (Cp) = " + u2_prime(1));
disp("u2 at x=0 (Cp) = " + u2(1) * 1000 + " mm");
disp("u2_prime at x=beam_length (Cg) = " + u2_prime(end));
disp("u2 at x=beam_length (Cg) = " + u2(end) * 1000 + " mm"); % should be 0 anyways

% u2" at sigma11

u2_pp = M3 ./ H;        % curvature (kappa) = u2" [1/m]

% representative x2 for each ply: mid-plane
x2_mid = 0.5*(x2_bounds(1:end-1) + x2_bounds(2:end));  % 1x8

% sigma at top and bottom for each ply and x location
sigma11_top = zeros(n_layers, length(x1_vec));        % [MPa]
sigma11_bottom = zeros(n_layers, length(x1_vec));     % [MPa]

for k = 1:n_layers
    Ex_MPa = Ex1x1_layers(k);

    % top fiber (theta = +90 deg): x2 = +r
    sigma11_top(k,:) = Ex_MPa * (u1prime_axial - (x2_mid(k)) * u2_pp);

    % bottom fiber (theta = -90 deg): x2 = -r
    sigma11_bottom(k,:) = Ex_MPa * (u1prime_axial + (x2_mid(k)) * u2_pp);
end

% Find global max absolute stress
absmax = 0;
xcrit_index = 1;
kcrit = 1;
sidecrit = "top";

for k = 1:n_layers
    [m1, i1] = max(abs(sigma11_top(k,:)));
    if m1 > absmax
        absmax = m1; xcrit_index = i1; kcrit = k; sidecrit = "top";
    end

    [m2, i2] = max(abs(sigma11_bottom(k,:)));
    if m2 > absmax
        absmax = m2; xcrit_index = i2; kcrit = k; sidecrit = "bottom";
    end
end

disp("Max |sigma11| = " + absmax + " MPa");
disp("Occurs at x1 = " + x1_vec(xcrit_index) + "m, ply k = " + kcrit + ", side = " + sidecrit);

% ================= KEVIN'S COMMENTS ================= %

% Governing Equation 1:
%
% Governing Equation 2:
%
% Boundary Conditions: (Investigating section of rocket from Cp -> Cg -> beam_length)
% @ Left Edge (Origin = O = Center of Pressure)
%
% @ Right Edge (x1 = beam_length)
%
% u2(R) = 0, u2'(R) = 0
% u2(L) = u2(R) @ x1=cpcg
% u2'(L) = u2'(R) @ x1=cpcg where u2' is the derivative of u2 wrt x1 (curvature)
%
% Weight force in middle of beam section causes discontinuity in shear and bending moment
% -> Must analyze the beam in two sections (x1=0->cpcg (left) and x1=cpcg->beam_length (right))
% Find the axial displacement u1 due to bending for both left and right sides
% Then find the left and right side total axial displacement by adding u1_axial to both u1(L/R)
% Find the strain and the stress associated with the axial displacement utot_L and utot_R
% Pick the higher stress value between L and R for further analysis
% -> Step by step layer analysis

% ==================================================== %
