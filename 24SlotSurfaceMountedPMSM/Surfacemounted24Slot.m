clear
close all
clc
%% This files is the analytical calculation for the 24Slot surfacemounted PMSM
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
d_st = 0.0008;              % Diamter of the strand [m]
L_stack = 0.150;            % Active lenght of the machine [m]
OD_st = 0.250;              % Stator outer diameter [m]
ID_st = 0.150;              % Stator inner diamter [m], Assumption
l_gap = 0.0008;             % Airgap [m]
t_PM = 0.005;               % Magnet thickness [m]
D_rot = ID_st - 2 * (l_gap + t_PM);     % Rotor diamter [m]

% Reluctance calculations

%PM
w_PM = (D_rot * pi * 0.995) / 8;        % Width of the PM [m]
A_PM = L_stack * w_PM;                  % Area of PM [m2]
R_PM = t_PM / (mu_0 * mu_PM * A_PM);    % Reluctance of the PM


% Air gap
w_gap = ((ID_st - l_gap) * pi) / 8; % Width of the airgag [m], the circumfrance of the 
A_gap = L_stack * w_gap;            % Area of the airgap [m2]
R_gap = l_gap  / (mu_0 * A_gap);    % Reluctance of the arigap


% Calculating the MMF of the PM
H_c = 920000;                       % Coercivity of the PM, choosing N36Z_20
MMF = H_c * t_PM;                   % MMF of the PM [H/m]

% Short circuit flux of PM
Phi_r = MMF / R_PM;                 % [Wb]

R_tot = 2 * (R_PM + R_gap);         % Total reluctance of the machine
MMF_tot = 2 * MMF;

% Calculating the flux
phi_noLoad = MMF_tot / R_tot;         % Short circuit flux [Wb], Assuming ideal steel with no MMF loss and infinit permeability 

% Calculating slot area
Hs0 = 0.004;
A_outer = (((2/3) * (OD_st - ID_st)) + ID_st)^2 * pi / 4;   % Assuming the slot covers 2/3 of the yoke area
A_iner = (ID_st + Hs0)^2 * pi / 4;
A_slot = (A_outer - A_iner) / (2 * Q);                      % Slot area [m2]
A_slot =  403.05*10^-6;                                                   

% Rated current
J = 5.3*10^6;                           % Curent density per strand [A/mm2]
N_series = 1;                           % Number of series branches
N_par = (2 * p) / N_series;             % Number parallel branches (N_par * N_ser = 2*p)
A_strand = (d_st^2 * pi) / 4;           % Strand area, neglecting the insulatioin

% Calculating no-load voltage
RPM = 1000;
f = (RPM  * p) / 60;   % Electrical rotation per second (rps)
fprintf("Calculating no-load voltage at: %d \n",RPM)

% Calculating number of turns and current as function of number of strands
n = 13;     % number of loops/ number of possible strands
N_st = linspace(1,n,n);
for i = 1:length(N_st)
    N_turn(i) = (A_slot * 0.45) ./ (A_strand * N_st(i));
    I_rms(i) = N_par * N_st(i) * A_strand * J;
    E_0(i) = (sqrt(12) * pi * f .* N_turn(i) * k_w * q * r * p * phi_noLoad) / N_par;
    fprintf('N_st = %d, N_turn = %d, I_phase = %.1f, E_0 = %0.1f \n', N_st(i), floor(N_turn(i)) ,I_rms(i),E_0(i));
end

% Choosing
N_turn = 60;
N_st = 5;
fprintf("Choosing the number of turns and strands:\n N_turn = %.0f, N_str = %.0f \n", N_turn, N_st);

% Magnetic flux linkage
Psi_PM = k_w * p * q * r * phi_noLoad * N_turn / N_par;

% d and q-axis reluctances
R_d = R_tot;   % d-axis reluctance
R_q = R_tot;   % q-axis reluctance
% Assuming the steel has no loss and is ideal. The reluctance is only due to the magnets and the air.
% since the machine is salient, the d and q-axis reluctances are the same.
N_d = N_turn * k_w * q * r;
N_q = N_d;
L_d = (N_d^2 * p) / (R_d * N_par^2);    % d-axis inductance
L_q = (N_q^2 * p) / (R_q * N_par^2);    % q-axis inductance


EMF = (sqrt(2) * pi * f .* N_turn * k_w * q * r * p * phi_noLoad) / N_par;     % RMS phase voltage
EMF_phase2phase_max = EMF * sqrt(6);    % maximum phase2phase voltage
I_rms_phase = N_par * N_st * A_strand * J;
fprintf('The maximume phase2phase induced voltage at: \n %.0f [RPM] is  %.2f [V] and I_phase = %.2f [A] \n', RPM, EMF_phase2phase_max, I_rms_phase)


%% Flux from simulation for 24 slot
DB24 = readtable('withAdot_bitch.csv');
% Extract columns as arrays (vectors)
D24 = DB24.('Distance_mm_');  % Adjust this to the actual column name from your CSV file
B24 = DB24.('B_normal__');
% Save the variables to a .mat file
save('vectors_for_matlab.mat', 'D24', 'B24');

delta = 0.5845 * 10^-3; % Distance step delta
phi_sim = sum(delta * L_stack .* -B24);  % [Wb]
Psi_PMSim = k_w * p * q * r * phi_sim * N_turn / N_par;

%% Calculating Ld and Lq inductances exported from simulation

LdLq = readtable('24slot_ldlq_updated.csv');
time = LdLq.("Time_ms_");
LdSim = mean(LdLq.("L_d_Axis__uH_"))*10^-6;
LqSim = mean(LdLq.("L_q_Axis__uH_"))*10^-6;

%% Flux linkages of the phases
psi_phases = readtable('24slot_fluxLinkage_updated.csv');
psi_time = psi_phases.Time_ms_;
psi_A = psi_phases.FluxLinkage_PhaseA__Wb_;
psi_B = psi_phases.FluxLinkage_PhaseB__Wb_;
psi_C = psi_phases.FluxLinkage_PhaseC__Wb_;
plot(psi_time,psi_A,'r'), hold on
plot(psi_time,psi_B,'b'), hold on
plot(psi_time,psi_C,'k'), grid on
legend('\Psi_A','\Psi_B','\Psi_C')
xlabel('time[ms]'), ylabel('\Psi [Wb]')

% %% Calculating L_d and L_q from simulations 48 slot
% wt = 0;
% T = (2/3)*[cos(wt) cos(wt-(2*pi/3)) cos(wt+(2*pi)/3)
% -sin(wt) -sin(wt-(2*pi)/3) -sin(wt+(2*pi)/3)
% 1/2 1/2 1/2];
% 
% L_phase = [53.674 -7.0635 -7.0335
%            -7.0635 53.28 -7.0315
%            -7.0335 -7.0315 53.519];
% L_dq = T * L_phase * inv(T);

