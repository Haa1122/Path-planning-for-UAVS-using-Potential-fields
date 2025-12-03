%% buildPotentialField.m
clear; clc; close all;

%% === 1. LOAD RISK MAP ===
risk = readmatrix('riskMap.csv');
lat  = readmatrix('lat.csv');
lon  = readmatrix('lon.csv');

[nRows, nCols] = size(risk);

%% === 2. USER-DEFINED START AND GOAL (pick any) ===
% Example: choose two points manually (will improve later)
start_i = 400; start_j = 200;   % row, col
goal_i  = 150; goal_j  = 700;

fprintf('Start = (%d,%d), Goal = (%d,%d)\n', start_i, start_j, goal_i, goal_j);

%% === 3. ATTRACTIVE POTENTIAL FIELD ===
k_att = 1e-3;  % small factor to avoid dominating

% Create meshgrid in index space
[jGrid, iGrid] = meshgrid(1:nCols, 1:nRows);

U_att = k_att * ((iGrid - goal_i).^2 + (jGrid - goal_j).^2);

%% === 4. REPULSIVE POTENTIAL FIELD FROM RISK MAP ===
% Normalize risk to [0,1] so values don't explode
risk_norm = risk / max(risk(:));

k_rep = 5;      % scaling factor (adjustable)
beta  = 2;      % risk exponent

% Repulsive potential
U_rep = k_rep * (risk_norm .^ beta);

%% === 5. TOTAL POTENTIAL ===
U = U_att + U_rep;

%% === 6. VISUALIZE THE POTENTIAL FIELD ===
figure;
imagesc(U);
axis equal tight;
set(gca,'YDir','normal');
colormap(jet); colorbar;
title('Total Potential Field U = U_{att} + U_{rep}');
xlabel('Column (j)');
ylabel('Row (i)');

hold on;
plot(start_j, start_i, 'wo', 'MarkerSize',10, 'LineWidth',2);
plot(goal_j, goal_i, 'wx', 'MarkerSize',12, 'LineWidth',3);
legend('','Start','Goal');
