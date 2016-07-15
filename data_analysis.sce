clear;
close;
close;
clc;

// Import Data
Tsim = fscanfMat('Times_SndSens.txt');
N_sim = length(Tsim);

T_aut = fscanfMat('Times_SndCom.txt');
N_aut = length(T_aut);

T_gs = fscanfMat('Times_GS.txt');
N_gs = length(T_gs);


// Compute the DT
DT_sim = Tsim(2 : N_sim) - Tsim(1 : N_sim - 1);
DT_aut = T_aut(2 : N_aut) - T_aut(1 : N_aut - 1);
DT_gs = T_gs(2 : N_gs) - T_gs(1 : N_gs - 1);

// Statistic
sdt_sim = stdev(DT_sim)
sdt_aut = stdev(DT_aut)
sdt_gs = stdev(DT_gs)

f1 = figure();
plot(DT_sim);

f2 = figure();
plot(DT_aut);

f3 = figure();
plot(DT_gs);
