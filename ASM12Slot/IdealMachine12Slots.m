%% ASM Task 12Slts idealized iron

clc
clear
clear all

%% Constants
mu_0 = 4 * pi * 10^-7;
%%
m = 3;
p = 2;
r = 1;                      % Number of winding layers
Q = 12;                     % Number of slots
q = Q / (2 * m * p);        % Slots per pole per pahse

%% Dimensioins
L_stack = 0.150;
R_gap = 0.0786;
L_gap = 0.0006;     % Airgap
W_tooth = 0.0364;   % Width of one tooth, assuming slot and tooth has the same width
W_gap2 = (2 * pi * R_gap) / (2 * p);

W_gap = 3 * W_tooth;
A_gap = W_gap2 * L_stack;


B_s0 = 0.00002;       % Slot gap
H_s0 = 0.0002;        % Width of the slot gap
%% Reluctance
R_gap = L_gap / (mu_0 * A_gap);
R_Leakage = 3 * (B_s0 / (mu_0 * H_s0 * L_stack));

%% Design parameters
 kw = 1;        % Assuming winding factor of 1
 N_turn = 20;       % Assuming
 N_ser = 1;
 N_par = (2 * p) / N_ser;       % Number of parallel branches
 N_eq = (kw * p * q * r * N_turn) / N_par;        % Equivalent number of turns
 
 %% Inductances

 L_m = (N_eq^2) / (R_gap);
 L_mLeakage =  N_eq^2 / R_Leakage;  
 
 InductanceRatio = L_mLeakage / L_m

 %% From simulations
 Psi_p2p0Slip = 3.594;       % Pk2Pk fluxlinkage of phase A
 Psi_0Slip = Psi_p2p0Slip / 2;
 T = 0.02;              % Period time
 Is_p2p0Slip = 707.107;      % Pk2Pk current of phase A
 Is_0Slip = Is_p2p0Slip/2;
 Es_p2p0Slip = 1124.732;     % Induced voltage
 Es_0Slip = Es_p2p0Slip / 2;

 f0 = 1 / T;
 w0 = 2 * pi * f0;
Ls_psi0Slip = Psi_0Slip / Is_0Slip;
% Ls_E0 = Es_0Slip / (w0 * Is_0Slip);

