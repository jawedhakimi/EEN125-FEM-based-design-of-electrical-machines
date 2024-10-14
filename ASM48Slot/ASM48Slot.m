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
p = 2;
q = Q / (2 * p * m);
k_fill = 0.4;

J = 30*10^6;        % Current density 30A/m^2
A_strand = pi * d_strand^2 / 4;
I_max = 250;        % 
%% 
N_ser = 1;
N_par = (2 * p) / N_ser;

N_strand = floor (I_max / (J * N_par * A_strand));
N_trun = floor ((k_fill * A_slot) / (A_strand * r * N_strand));

