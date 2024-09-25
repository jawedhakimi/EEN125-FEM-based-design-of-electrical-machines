
%% Part two of the task 48slot V-shape magnets
clear
close all
clc

% Constants 
mu_0 = 4* pi *10^-7; % Vacuum permeability [H/m]
mu_PM = 1.03;        % Magnets relative permeability


% Given parameters
Q = 48;             % Number of slots
p = 4;              % Number of pole pairs
k_w = 1;            % Winding factor, assuming 1
m = 3;              % Number of phases
r = 1;              % winding layer, 1
q = Q / (2 * p * m);              % Slots per pole per phaseygv8ytvygvyt

% Dimensions
d_st = 0.8;         % Diamter of the strand [mm]
l_stack = 150;      % Active lenght of the machine [mm]
OD_st = 269.24;        % Stator outer diameter [mm]
ID_st = 161.9;        % Stator inner diamter [mm], Assumption
OD_rot = 160.4;
l_gap = 0.73;        % Airgap [mm]


w_PM = 32;          % Magets width
t_PM = 6.48;        % Magnet thickness [mm]
D_rot = ID_st - 2 * l_gap;     % Rotor diamter [mm]


Hs0 = 1.03;
Hs1 = 0;
Hs2 = 29.5;
Rs = 5;
%% MMF of the PM
H_c = 920000;               % Coercivity of the PM, choosing N36Z_20
F_c = H_c * t_PM * 10^-3;   % F_c of the PM [H/m]

% Short circuit flux of PM
A_PM = l_stack * w_PM;
B_r = mu_0 * mu_PM * H_c;
Phi_r = B_r * A_PM * 10^-6;     % [Wb]


%% Reluctance calculations
%PM

A_PM = l_stack * w_PM;                        % Area of PM [mm2]
R_PM = (t_PM * 10^-3) / (mu_0 * mu_PM * A_PM * 10^-6);      % Reluctance of the PM
%%
% Stator
load BH_data_Prius_tut.mat
% Extract H and B vectors from the data
H_steel = BH_data(:, 1);  % First column corresponds to H (A_per_meter)
B_steel = BH_data(:, 2);  % Second column corresponds to B (tesla)
mu_steel = B_steel ./ (mu_0 .* H_steel);
%mu_steel = B_steel(9) / (mu_0 * H_steel(9));

% figure(2)
% plot( - H_steel, B_steel, 'r',LineWidth = 2)
% xlabel('H'), ylabel('B')
%%
% Reluctance of the tooth
w_t = 5.73;     % Width of the tooth
l_tooth = Hs0 + Hs1 + Hs2 + Rs;     % [mm] length of the flux path in the stator tooth
A_tooth = 2.5 * w_t * l_stack;      % [mm] width of the stator yoke
R_tooth = (l_tooth * 10^-3) ./ (mu_0 .* mu_steel .* A_tooth .* 10^-6);  % Reluctance of the tooth

% Reluctance of the stator yoke
dsy1 = ID_st + 2 * (Hs0 + Hs1 + Hs2 + Rs);      % Inner diameter of the flux path on the stator yoke
dsy2 = OD_st;       % outer diameter of the flux path in the stator yoke
dsym = (dsy2 + dsy1) / 2;                   % [mm] middle diamter of the stator yoke     
Asy = ((dsy2 - dsy1) / 2) * l_stack;        % [mm2]the crossectional area of the flux path in the stator yoke
lsy = ((4 * pi * dsym) / 48) + ((dsy2 - dsy1) / 2);     % [mm] length of the flux path in the stator, taking four teeth
R_st_yoke = lsy .* 10^-3 ./ (mu_0 .* mu_steel .* Asy .* 10^-6);    % reluctance of the stator yoke

% Air gap            % Width of the airgag [mm], the circumfrance of the
l_gap = 0.73;        % Airgap [mm]
A_gap = (pi * (ID_st + OD_st) * l_stack) / (16 * 2);        % Area of the airgap [mm2]
R_gap = (l_gap * 10 ^-3) / (mu_0 * A_gap * 10^-6);          % Reluctance of the arigap

% Rotor
l_rot = ((pi / 2) * (4 * pi * OD_rot) / 48) - 2 * t_PM;         % [mm] length of the flux path in the stator yoke
w_rot = w_PM;                                                   % [mm] width of the stator yoke
A_rot  = w_PM * l_stack;                                        % Crossection area of rotor
R_rot_yoke = l_rot .* 10^-3 ./ (mu_0 .* mu_steel .* w_rot .* l_stack .* 10^-6);    % reluctance of the stator yoke

% Total reluctance
R_tot = 2 * (R_tooth + R_gap + R_PM) + R_st_yoke + R_rot_yoke;
phi = F_c ./ R_tot;


%% Slot area
Hs0 = 4;
A_outer = (((2/3) * (OD_st - ID_st)) + ID_st)^2 * pi / 4;   % Assuming the slot covers 2/3 of the yoke area and the slot has the same area as the teeth
A_iner = (ID_st + Hs0)^2 * pi / 4;
A_slot = (A_outer - A_iner) / (2 * Q);                           % Slot area [mm2]

% Rated current
J = 5;                          % Curent density per strand [A/mm2]
N_par = 8;                      % Number parallel branches
A_strand = 0.8^2 * pi / 4;      % Strand area, neglecting the insulatioin
