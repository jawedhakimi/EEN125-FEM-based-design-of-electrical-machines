
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
m = 10^-3;
d_st = 0.8 * m;         % Diamter of the strand [m]
L_stack = 150 * m;      % Active lenght of the machine [m]
OD_st = 269.24 * m;        % Stator outer diameter [m]
ID_st = 161.9 * m;        % Stator inner diamter [m]
OD_rot = 160.4 * m;
l_gap = 0.73 * m;        % Airgap [m]


w_PM = 32 * m;          % Magets width[m]
t_PM = 6.48 * m;        % Magnet thickness [m]
D_rot = (ID_st - 2 * l_gap) * m;     % Rotor diamter [m]


Hs0 = 1.03 * m;
Hs1 = 0 * m;
Hs2 = 29.5 * m;
Rs = 5 * m;
%% MMF of the PM
H_c = 920000;               % Coercivity of the PM, choosing N36Z_20
F_c = H_c * t_PM;   % F_c of the PM [H/m]

% Short circuit flux of PM
A_PM = L_stack * w_PM;
B_r = mu_0 * mu_PM * H_c;
Phi_r = B_r * A_PM;     % [Wb]


%% Reluctance calculations
%PM

A_PM = L_stack * w_PM;                    % Area of PM [m2]
R_PM = t_PM / (mu_0 * mu_PM * A_PM);      % Reluctance of the PM

% Stator
load BH_data_Prius_tut.mat
% Extract H and B vectors from the data
H_steel = BH_data(:, 1);  % First column corresponds to H (A_per_meter)
B_steel = BH_data(:, 2);  % Second column corresponds to B (tesla)
mu_steel = B_steel ./ (H_steel .* mu_0);

% Reluctance of the tooth
w_t = 5.73 * m;     % Width of the tooth
l_tooth = Hs0 + Hs1 + Hs2 + Rs;     % [m] length of the flux path in the stator tooth
A_tooth = 2.5 * w_t * L_stack;      % [m] width of the stator yoke
R_tooth = l_tooth ./ (mu_0 .* mu_steel .* A_tooth);  % Reluctance of the tooth

% Reluctance of the stator yoke
dsy1 = ID_st + 2 * (Hs0 + Hs1 + Hs2 + Rs);      % Inner diameter of the flux path on the stator yoke
dsy2 = OD_st;       % outer diameter of the flux path in the stator yoke
dsym = (dsy2 + dsy1) / 2;                   % [m] middle diamter of the stator yoke     
Asy = ((dsy2 - dsy1) / 2) * L_stack;        % [m2]the crossectional area of the flux path in the stator yoke
lsy = ((4 * pi * dsym) / 48) + ((dsy2 - dsy1) / 2);     % [m] length of the flux path in the stator, taking four teeth
R_st_yoke = lsy ./ (mu_0 .* mu_steel .* Asy);    % reluctance of the stator yoke

% Air gap
A_gap = (pi * (ID_st + OD_st) * L_stack) / (16 * 2);        % Area of the airgap [m2]
R_gap = l_gap / (mu_0 * A_gap);          % Reluctance of the arigap

% Rotor
l_rot = ((pi / 2) * (4 * pi * OD_rot) / 48) - 2 * t_PM;         % [m] length of the flux path in the stator yoke
w_rot = w_PM;                                                   % [m] width of the stator yoke
A_rot  = w_PM * L_stack;                                        % Crossection area of rotor
R_rot_yoke = l_rot ./ (mu_0 .* mu_steel .* w_rot .* L_stack);    % reluctance of the stator yoke

% Total reluctance
R_tot = 2 .* (R_tooth + R_gap + R_PM + R_rot_yoke) + R_st_yoke ;
phi_load = linspace (0,Phi_r,25);
F_load =  (phi_load .* transpose(R_tot))';


%%
%PLOT OF LOAD & SOURCE LINE
%assuming B_air gap
B_gap = 0.01:0.001:2;
phi_gap = B_gap * A_gap;
%we can compute flux densities in different parts of the iron core
Brot=phi_gap/A_rot;
Bstator_yoke=phi_gap/Asy;
Bstator_tooth=phi_gap/A_tooth;
%now we have to compute the H values in order then to get mmfs
%interpolation must be used since probably u do not get an exact value :(
method = 'linear';
H_gap = B_gap/mu_0;
Hstator_tooth=interp1(B_steel,H_steel,Bstator_tooth,method);
Hstator_yoke=interp1(B_steel,H_steel,Bstator_yoke,method);
Hrot=interp1(B_steel,H_steel,Brot,method);

%Let's start with mmf drops
Fr=Hrot*l_rot;
Fsy=Hstator_yoke*lsy;
Fst=Hstator_tooth*l_tooth;
Fag=H_gap*l_gap;

MMFload=Fr+2*Fag+2*Fst+Fsy;

%source line
Bmag=[0 B_r];
Hmag=[-H_c 0];
MMFmag=2*Hmag*t_PM;
Phimag=Bmag*A_PM;

figure(2)
plot(MMFmag,Phimag,'r','LineWidth',2)
xlim([min(MMFmag) 0]),hold on
plot(-MMFload,phi_gap,'m','LineWidth',2), grid on
xlabel('MMF'), ylabel('Flux')
%defining some boundaries to better understanding
xlim([min(MMFmag) 0]);
ylim([0 Phi_r]);

%% Choose operating points flux from the fgure
phi = 0.0040888;      % [Wb]

% Slot area
A_slot = (336.8 + 66.24) * 10^-6;                           % Slot area [m2]

% Rated current
J = 5 * 10^6;                          % Curent density per strand [A/m2]
N_par = 8;                      % Number parallel branches
A_strand = d_st^2 * pi / 4;      % Strand area, neglecting the insulatioin

% Calculating number of turns and current as function of number of strands
n = 12;     % number of loops/ number of possible strands
N_st = linspace(1,n,n);
for i = 1:length(N_st)
    N_turn(i) = (A_slot * 0.45) ./ ((pi / 4) * d_st^2 * N_st(i));
    I_rms(i) = N_par * N_st(i) * A_strand * J;
    fprintf('N_st = %.0f, N_turn = %.0f, I_phase = %.1f,\n', N_st(i), floor(N_turn(i)) ,I_rms(i));
    i = i+1;
end


fprintf("Calculating no-load voltage \n")
% Calculating no-load voltage
RPM_max = 14000;    % amximum speed of the machine
RPM_rated = RPM_max * (1 / 3);  % Rated speed of the machine, just before field weakening


RPM = 5500;
f = (RPM  * p) / 60;   % Electrical rotation per second (rps)

fprintf("The speed is: = %.0f [rpm] \n", RPM)
for i = 1:length(N_turn)
    E_0(i) = (sqrt(12) * pi * f .* N_turn(i) * k_w * q * r * p * phi) / N_par;
    fprintf("Induced phase2phase voltage = %.1f, N_turn = %.1f\n", E_0 (i), floor(N_turn(i)));
    i = i+1;
end

% Choosing 
N_turn = 68;    % Chose from the table
N_st = 5;
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


EMF = (sqrt(2) * pi * f .* N_turn * k_w * q * r * p * phi) / N_par;     % RMS phase voltage
EMF_phase2phase_max = EMF * sqrt(6);    % maximum phase2phase voltage
I_rms_phase = N_par * N_st * A_strand * J;
fprintf('The maximume phase2phase induced voltage at: \n %.0f [RPM] is  %.2f [V] and I_phase = %.2f [A] \n', RPM, EMF_phase2phase_max, I_rms_phase)


