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
L_gap = 0.0005;     % Airgap
W_tooth = 0.0196;   % Width of one tooth, assuming slot and tooth has the same width
W_gap = 3 * W_tooth;
A_gap = W_gap * L_stack;

B_s0 = 0.002;       % Slot gap
H_s0 = 0.002;       % Width of the slot gap
%% Reluctance
R_gap = L_gap / (mu_0 * A_gap);
R_Leakage = 3 * (B_s0 / (mu_0 * H_s0 * L_stack)) + 2 * R_gap;

%% Design parameters
 kw = 1;        % Assuming winding factor of 1
 N_turn = 20;       % Assuming
 N_ser = 1;
 N_par = (2 * p) / N_ser;       % Number of parallel branches
 N_eq = kw * p * q * r * N_turn / N_par;        % Equivalent number of turns
 
 %% Inductances

 L_m = N_eq^2 / (2 * R_gap);
 L_mLeakage =  N_eq^2 / R_Leakage;  
 
 InductanceRatio = L_mLeakage / L_m