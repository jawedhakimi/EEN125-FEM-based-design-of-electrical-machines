

clear
close all
clc

% Constants 
mu_0 = 4* pi *10^-7; % Vacuum permeability [H/m]
mu_PM = 1.03;        % Magnets relative permeability


% Given parameters
Q = 24;             % Number of slots
p = 4;              % Number of pole pairs
k_w = 1;            % Winding factor, assuming 1
m = 3;              % Number of phases
r = 1;              % winding layer, 1
q = Q / (2 * p * m);              % Slots per pole per phase

% Dimensions
d_st = 0.8;         % Diamter of the strand [mm]
l = 150;            % Active lenght of the machine [mm]
OD_st = 250;        % Stator outer diameter [mm]
ID_st = 150;        % Stator inner diamter [mm], Assumption
l_gap = 0.8;        % Airgap [mm]
t_PM = 5;           % Magnet thickness [mm]
D_rot = ID_st - 2 * (l_gap + t_PM);     % Rotor diamter [mm]


% Reluctance calculations
w_PM = (D_rot * pi * 0.995) / 8;        % Width of the PM [mm]
A_PM = l * w_PM;                        % Area of PM [mm2]
R_PM = (t_PM * 10^-3) / (mu_0 * mu_PM * A_PM * 10^-6);      % Reluctance of the PM



w_gap = ((ID_st - l_gap) * pi) / 16;                         % Width of the airgap [mm]
A_gap = l * w_gap;                                          % Area of the airgap [mm2]
R_gap = (l_gap * 10 ^-3) / (mu_0 * A_gap * 10^-6);          % Reluctance of the arigap


% Calculating the MMF of the PM
H_c = 920000;               % Coercivity of the PM, choosing N36Z_20
MMF = H_c * t_PM * 10^-3;   % MMF of the PM [H/m]

% Short circuit flux of PM
Phi_r = MMF / R_PM;     % [Wb]

% Calculating the flux
phi = MMF / (R_PM + R_gap);         % Short circuit flux [Wb], Assuming ideal steel with no MMF loss and infinit permeability 

% Calculating slot area
A_outer = (((2/3) * (OD_st - ID_st)) + ID_st)^2 * pi / 4;   % Assuming the slot covers 2/3 of the yoke area
A_iner = 154^2 * pi / 4;
A_slot = (A_outer - A_iner) / 48;                           % Slot area [mm2]

% Rated current
J = 5;                          % Curent density per strand [A/mm2]
N_par = 8;                      % Number parallel branches
A_strand = 0.8^2 * pi / 4;      % Strand area, neglecting the insulatioin

% Calculating number of turns and current as function of number of strands
n = 13;     % number of loops/ number of possible strands
N_st = linspace(1,n,n);
for i = 1:length(N_st)
    N_turn(i) = (A_slot * 0.45) ./ ((pi / 4) * d_st^2 * N_st(i));
    I_rms(i) = N_par * N_st(i) * A_strand * J;
    %fprintf('N_st(%d) = %.0f, N_turn(%d) = %.0f, I_phase(%d) = %.1f,\n', i, N_st(i), i, floor(N_turn(i)), i ,I_rms(i));
    i = i+1;
end


fprintf("Calculating no-load voltage \n")



% Calculating no-load voltage
RPM_max = 14000;    % amximum speed of the machine
RPM_rated = RPM_max * (1 / 3);  % Rated speed of the machine, just before field weakening


RPM = 5000;
f = (RPM  * p) / 60;   % Electrical rotation per second (rps)

fprintf("The speed is: = %.0f [rpm] \n", RPM)
for i = 1:length(N_turn)
    E_0(i) = (sqrt(6) * pi * f .* N_turn(i) * k_w * q * r * p * phi) / N_par;
    %fprintf("Induced phase2phase voltage (%d) = %.1f, N_turn(%d) = %.1f\n", i, E_0 (i), i, floor(N_turn(i)));
    i = i+1;
end

% Choosing 
N_turn = 56;    % Chose from the table
N_st = 6;
fprintf("Choosing the number of turns and strands:\n N_turn = %.0f, N_str = %.0f \n", N_turn, N_st);

% Magnetic flux linkage
Psi_PM = k_w * p * q * r * phi * N_turn / N_par;

% d and q-axis reluctances
R_d = (R_PM + R_gap);   % d-axis reluctance
R_q = (R_PM + R_gap);   % q-axis reluctance
% Assuming the steel has no loss and is ideal. The reluctance is only due to the magnets and the air.
% since the machine is salient, the d and q-axis reluctances are the same.
N_d = N_turn * k_w * q * r;
N_q = N_d;
L_d = (N_d^2 * p) / (R_d * N_par^2);    % d-axis inductance
L_q = (N_q^2 * p) / (R_q * N_par^2);    % q-axis inductance

% w = 2 * pi * (RPM * p / 60);    % [rad/s]
% U_lineAmpMax = 800 * 0.995;     % Terminal voltage [V]
% N_turnNew = (U_lineAmpMax * N_par) / (sqrt(3) * k_w * p * q * phi * w * r * p)


EMF = (sqrt(2) * pi * f .* N_turn * k_w * q * r * p * phi) / N_par;     % RMS phase voltage
EMF_phase2phase_max = EMF * sqrt(6);    % maximum phase2phase voltage
I_rms_phase = N_par * N_st * A_strand * J;
fprintf('The maximume phase2phase induced voltage at: \n %.0f [RPM] is  %.2f [V] and I_phase = %.2f [A] \n', RPM, EMF_phase2phase_max, I_rms_phase)

