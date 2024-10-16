clc
clear 

%% Constants
mu_0 = 4 * pi * 10^-7;

%% Dimensions
d_strand = 0.0008; % Strand diameter [m]
A_slot = 202 * 10^-6;


%% Electrical parameters
Q = 48;
m = 3;
r = 1;
N_pp = 2;
N_p = 4;
q = Q / (2 * N_pp * m);
k_fill = 0.4;
T_cu = 20;
Rho_Co = 16.8*10^-9*(1+3.9*10^-3 *(T_cu - 20 ));
N_turn = 17;
L_stack = 0.15;


J = 30*10^6;        % Current density 30A/m^2
A_strand = pi * d_strand^2 / 4;
I_max = 250;        % 
%% 
N_ser = 1;
N_pb = (2 * N_pp) / N_ser;

N_strand = floor (I_max / (J * N_pb * A_strand));
N_trun = floor ((k_fill * A_slot) / (A_strand * r * N_strand));


%% From Simulation

% Noload
f0_NoLoad = 50; % Hz
PsiA_NoLoadAbs = 0.7176; % Wb
iA_NoLoadAbs = 354;
uA_NoLoadAbs = 218;
Ls_NoLoad = PsiA_NoLoadAbs / iA_NoLoadAbs;

% Locked rotor
f0_LockedRotor = 12.5; %Hz
uA_LockedRotorAbs = 7.3852;
ThetaUA_LockedRotor = deg2rad(47.07);    % Deg
ThetaIA_LockedRotor = 0;        % Deg
iA_LockedRotorAbs = 354;

Ls_lambda = (uA_LockedRotorAbs * sin(ThetaUA_LockedRotor)) / (4 * pi * f0_LockedRotor * iA_LockedRotorAbs);
Lm = Ls_NoLoad - Ls_lambda;
Rr = (uA_LockedRotorAbs * cos(ThetaUA_LockedRotor)) / iA_LockedRotorAbs;

%% Stator resistance
% Rs = (Rho_Co * N_p * N_turn * q * L_stack) / (2 * A_strand * N_pb^2);
Np = 4;
m = 3;
Nlayer = 1;
Nturn = 17;
q = Q / (2 * Np * m);
Npb = 4;
Rslotcenter = 0.08693;
Npitch = 10;
Lpitch = (1.2 * 2 * pi * Rslotcenter * Npitch) / Q;
Lendext = 0.025;
Lew = 2 * (2 * Lendext + Lpitch);
Law = 2 * L_stack;
Lcoil = Law + Lew;
Astrand = (d_strand / 2)^2 * pi;
Acoil = Astrand * N_strand;
Rs = (Np * Nlayer * q * Rho_Co * Nturn^2 * Lcoil) / (2 * N_pb^2 * Acoil);

%% 

Rs2 = (3110 / 250^2) / 3;



