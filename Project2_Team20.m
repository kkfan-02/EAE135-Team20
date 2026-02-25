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
fprintf('H33_C = %.4f MPa*m^4  (= %.4e N*m^2)\n', H33_C, H33_C * 1e6);


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

disp("S_axial = " + S_axial + " N");
disp("u1prime_axial = " + u1prime_axial);
disp("u1_axial at x=0 (fixed end) = " + u1_axial(1) * 1000 + " mm");

% ***************** BENDING PROBLEM *****************
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

% Lift and weight forces in x2 direction
% Lift acts at x1=0 (center of pressure, free end)
% Weight acts at x1=cpcg (center of gravity)
% Beam is clamped at x1=beam_length
lift_x2   = lift * cos(deg2rad(AoA));          % [N] transverse component of lift
weight_x2 = Weight * g * cos(deg2rad(AoA));    % [N] transverse component of weight force

% -------------------------------------------------------------------------
% BENDING MOMENT M3(x1)
% -------------------------------------------------------------------------
% Sign convention: positive M3 causes positive curvature (u2'' = -M3/H33_C).
% Taking moments about the cut, looking LEFT, with upward forces positive:
%
%   Section L (0 <= x1 <= cpcg):
%     Only lift_x2 acts to the left of the cut.
%     M3_L(x1) = lift_x2 * x1
%
%   Section R (cpcg <= x1 <= beam_length):
%     Both lift_x2 and weight_x2 act to the left of the cut.
%     M3_R(x1) = lift_x2 * x1 - weight_x2 * (x1 - cpcg)
%
% Note: weight_x2 opposes lift_x2 (acts downward), hence the minus sign.

% Section L: x1 from 0 to cpcg
x1_L = linspace(0, cpcg, 200);
M3_L = lift_x2 * x1_L;                                    % [N*m]

% Section R: x1 from cpcg to beam_length
x1_R = linspace(cpcg, beam_length, 400);
M3_R = lift_x2 * x1_R - weight_x2 * (x1_R - cpcg);       % [N*m]

% Combined x1 and M3 vectors for plotting
x1_all = [x1_L, x1_R(2:end)];
M3_all = [M3_L, M3_R(2:end)];

fprintf('Max |M3| in Section L = %.4e N*m  at x1 = %.4f m\n', max(abs(M3_L)), x1_L(end));
fprintf('Max |M3| in Section R = %.4e N*m\n', max(abs(M3_R)));
[M3_max_val, idx_max] = max(abs(M3_all));
x1_at_M3max = x1_all(idx_max);
fprintf('Overall max |M3| = %.4e N*m  at x1 = %.4f m\n', M3_max_val, x1_at_M3max);

% -------------------------------------------------------------------------
% BENDING DISPLACEMENT u2(x1)
% -------------------------------------------------------------------------
% Governing ODE:   H33_C * u2''(x1) = -M3(x1)
%   => u2''(x1) = -M3(x1) / H33_C_SI
%
% Boundary conditions at the clamp (x1 = beam_length):
%   u2(beam_length)  = 0
%   u2'(beam_length) = 0
%
% Integrate from right to left using cumulative trapezoid so BCs are
% naturally satisfied at x1 = beam_length.

H33_C_SI = H33_C * 1e6;   % convert MPa*m^4 -> N*m^2

% --- Section R (integrate right-to-left, flip then flip back) ---
% Flip so integration runs from beam_length down to cpcg
x1_R_flip  = fliplr(x1_R);
M3_R_flip  = fliplr(M3_R);
u2pp_R_flip = -M3_R_flip / H33_C_SI;   % u2''

% First integration: u2'  (BC: u2'(beam_length)=0, i.e. starting value=0)
% Negate because cumtrapz with decreasing x returns -(true integral)
u2p_R_flip = -cumtrapz(x1_R_flip, u2pp_R_flip);   % running from beam_length leftward

% Second integration: u2  (BC: u2(beam_length)=0, i.e. starting value=0)
u2_R_flip  = -cumtrapz(x1_R_flip, u2p_R_flip);

% Flip back to left-to-right order
u2p_R = fliplr(u2p_R_flip);
u2_R  = fliplr(u2_R_flip);

% u2'(cpcg) and u2(cpcg) from Section R serve as BCs for Section L
u2p_R_at_cpcg = u2p_R(1);   % u2'(cpcg) from right section
u2_R_at_cpcg  = u2_R(1);    % u2(cpcg)  from right section

% --- Section L (integrate right-to-left, with matching BCs at x1=cpcg) ---
x1_L_flip  = fliplr(x1_L);
M3_L_flip  = fliplr(M3_L);
u2pp_L_flip = -M3_L_flip / H33_C_SI;

% First integration with BC u2'(cpcg) = u2p_R_at_cpcg
u2p_L_flip = u2p_R_at_cpcg - cumtrapz(x1_L_flip, u2pp_L_flip);

% Second integration with BC u2(cpcg) = u2_R_at_cpcg
u2_L_flip  = u2_R_at_cpcg  - cumtrapz(x1_L_flip, u2p_L_flip);

u2p_L = fliplr(u2p_L_flip);
u2_L  = fliplr(u2_L_flip);

% Combined displacement
u2_all  = [u2_L,  u2_R(2:end)];
u2p_all = [u2p_L, u2p_R(2:end)];

fprintf('Max transverse displacement u2 = %.4f mm  at x1 = 0 (free end)\n', u2_L(1)*1000);

% =========================================================================
% AXIAL DISPLACEMENT FROM BENDING  u1_bend(x1) = -x2 * u2'(x1)
% (evaluated at the outermost fiber, x2 = r_outer, worst case)
% =========================================================================
r_outer = x2_bounds(end);   % outermost radius [m]

u1_bend_all = -r_outer * u2p_all;   % [m], worst-case fiber

% Total axial displacement (axial + bending contribution at outermost fiber)
% Interpolate u1_axial onto the combined x1 grid
u1_axial_all = u1prime_axial * (x1_all - beam_length);  % [m]
u1_total_all = u1_axial_all + u1_bend_all;               % [m]

% =========================================================================
% AXIAL STRESS  sigma1(x1, x2) LAYER-BY-LAYER
% =========================================================================
% From the PDF (page 6):
%   sigma1(x1, x2) = -(Ex1x1_k * x2 * M3(x1)) / H33_C
%
% Units: [MPa] * [m] * [N*m] / [MPa*m^4] = [N/m^2 * m * N*m / (N/m^2 * m^4)]
%      = [MPa*m*N*m] / [MPa*m^4] = [N*m^2/m^3] ... careful:
%
% Keeping H33_C in MPa*m^4 and M3 in N*m:
%   sigma1 [MPa] = -(Ex_MPa [MPa] * x2 [m] * M3 [N*m]) / (H33_C [MPa*m^4])
%                = -(Ex * x2 * M3) / H33_C   [MPa * m * N*m / MPa*m^4]
%                = -(Ex * x2 * M3) / H33_C   [N*m^2 / m^4] = [N/m^2] = [Pa]
%
% To get MPa: divide M3 by 1e6 (convert N*m -> MN*m = MPa*m^3... cleaner below):
%   sigma1 [MPa] = -(Ex_MPa * x2 * M3_Nm) / (H33_C_MPa_m4 * 1e6)
%
% Simplest: work in SI throughout for stress, then convert result to MPa.

% Find M3 at location of maximum |M3| (already found above)
M3_for_stress = M3_all(idx_max);   % [N*m]  (signed)

fprintf('\n--- Axial Stress at x1 = %.4f m (max |M3| location) ---\n', x1_at_M3max);
fprintf('M3 = %.4e N*m\n', M3_for_stress);
fprintf('H33_C = %.4e MPa*m^4\n', H33_C);
fprintf('%-6s  %-8s  %-12s  %-12s  %-12s  %-14s  %-14s  %-14s\n', ...
    'Layer', 'Angle', 'x2_top[m]', 'sig_top[MPa]', 'sig_bot[MPa]', 'Tensile Allow', 'Compressive Allow', 'Status');

% Allowables table (tensile positive, compressive negative) [MPa]
allowables_tensile     = [StrengthTensile_0deg,     StrengthTensile_45deg,  ...
                          StrengthTensile_45deg,     StrengthTensile_90deg,  ...
                          StrengthTensile_90deg,     StrengthTensile_45deg,  ...
                          StrengthTensile_45deg,     StrengthTensile_0deg];
allowables_compressive = [StrengthCompressive_0deg,  StrengthCompressive_45deg, ...
                          StrengthCompressive_45deg,  StrengthCompressive_90deg, ...
                          StrengthCompressive_90deg,  StrengthCompressive_45deg, ...
                          StrengthCompressive_45deg,  StrengthCompressive_0deg];

sigma_top_upper = zeros(1, n_layers);   % stress at top face, upper half
sigma_bot_upper = zeros(1, n_layers);   % stress at bottom face, upper half

first_fail_layer = 0;
first_fail_angle = 0;

for k = 1:n_layers
    x2_top = x2_bounds(k+1);   % outer face of layer k [m]
    x2_bot = x2_bounds(k);     % inner face of layer k [m]
    Ex_k   = Ex1x1_layers(k);  % [MPa]
    angle_k = ply_angles_deg(k);

    % Stress in MPa (units: MPa * m * N*m / (MPa * m^4) = N/m^2 = Pa -> /1e6 for MPa)
    sig_top = -(Ex_k * x2_top * M3_for_stress) / (H33_C * 1e6);  % [MPa]
    sig_bot = -(Ex_k * x2_bot * M3_for_stress) / (H33_C * 1e6);  % [MPa]

    sigma_top_upper(k) = sig_top;
    sigma_bot_upper(k) = sig_bot;

    % Apply factor of safety:
    %   For tensile stress: allowable = strength / FoS
    %   For compressive stress: allowable = strength / FoS  (strength already negative, so dividing makes it LESS negative = stricter)
    allow_T = allowables_tensile(k)     / FoS;
    allow_C = allowables_compressive(k) / FoS;   % divide by FoS to make it less negative (stricter)

    % Check top face
    if sig_top >= 0
        pass_top = sig_top <= allow_T;
        status_top = 'T';
    else
        pass_top = sig_top >= allow_C;
        status_top = 'C';
    end

    % Check bottom face (opposite sign due to symmetry)
    sig_bot_lower = -sig_top;   % lower half-plane mirror
    if sig_bot_lower >= 0
        pass_bot = sig_bot_lower <= allow_T;
    else
        pass_bot = sig_bot_lower >= allow_C;
    end

    pass_layer = pass_top && pass_bot;
    if pass_layer
        status_str = 'PASS';
    else
        status_str = '*** FAIL ***';
        if first_fail_layer == 0
            first_fail_layer = k;
            first_fail_angle = angle_k;
        end
    end

    fprintf('%-6d  %-8d  %-12.5f  %-12.4f  %-12.4f  %-14.2f  %-14.2f  %s\n', ...
        k, angle_k, x2_top, sig_top, sig_bot, allow_T, allow_C, status_str);
end

% Lower half-plane (mirrored): stresses are opposite sign
fprintf('\n--- Lower half-plane (x2 negative, mirrored stresses) ---\n');
fprintf('%-6s  %-8s  %-12s  %-12s\n', 'Layer', 'Angle', 'sig_top[MPa]', 'sig_bot[MPa]');
for k = 1:n_layers
    fprintf('%-6d  %-8d  %-12.4f  %-12.4f\n', ...
        k, ply_angles_deg(k), -sigma_top_upper(k), -sigma_bot_upper(k));
end

% =========================================================================
% SUMMARY
% =========================================================================
fprintf('\n========== SUMMARY ==========\n');
fprintf('H33_C             = %.4e N*m^2\n', H33_C_SI);
fprintf('S_axial           = %.4e N\n',     S_axial);
fprintf('Axial strain      = %.4e (dimensionless)\n', u1prime_axial);
fprintf('Max |M3|          = %.4e N*m  at x1 = %.4f m\n', M3_max_val, x1_at_M3max);
fprintf('Max u2 (tip)      = %.4f mm\n', u2_L(1)*1000);
if first_fail_layer > 0
    fprintf('First failing layer: %d  (angle = %d deg)\n', first_fail_layer, first_fail_angle);
else
    fprintf('All layers PASS under FoS = %.2f\n', FoS);
end

% =========================================================================
% PLOTS
% =========================================================================
figure(1);
plot(x1_all, M3_all/1e6, 'b-', 'LineWidth', 1.5);
xlabel('x1 [m]'); ylabel('M3 [MN*m]');
title('Bending Moment M3 along beam'); grid on;
xline(cpcg, '--r', 'CG'); xline(beam_length, '--k', 'Clamp');

figure(2);
plot(x1_all, u2_all*1000, 'r-', 'LineWidth', 1.5);
xlabel('x1 [m]'); ylabel('u2 [mm]');
title('Transverse Displacement u2 along beam'); grid on;
xline(cpcg, '--r', 'CG'); xline(beam_length, '--k', 'Clamp');
