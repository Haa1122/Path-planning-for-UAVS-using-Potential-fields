%% Gradient Descent Path Planning on Potential Field
clear; clc; close all;

%% === 1. Load Potential Field and Risk Map ===
risk = readmatrix('riskMap.csv');
lat  = readmatrix('lat.csv');
lon  = readmatrix('lon.csv');

% Same start and goal as before
start_i = 400; start_j = 200;
goal_i  = 150; goal_j  = 700;

[nRows, nCols] = size(risk);

% Create meshgrid for indices
[jGrid, iGrid] = meshgrid(1:nCols, 1:nRows);

%% === 3. ATTRACTIVE POTENTIAL FIELD (weaker) ===
k_att = 1e-5;
U_att = k_att * ((iGrid - goal_i).^2 + (jGrid - goal_j).^2);

%% === 4. REPULSIVE POTENTIAL FIELD USING LOG-RISK (stronger) ===
risk_clipped = max(risk, 1e-12);     % avoid log(0)
risk_log = -log(risk_clipped);       % convert risk to repulsion
risk_log = risk_log / max(risk_log(:));  % normalize

k_rep = 50;   % MUCH stronger repulsive field
U_rep = k_rep * risk_log;

%% === 5. TOTAL POTENTIAL ===
U = U_att + U_rep;
%% === 2. Gradient Descent Parameters ===
maxSteps = 5000;
stepSize = 1;       % move by one grid cell
tol = 2;            % stop if close to goal

path_i = start_i;
path_j = start_j;

path = [path_i, path_j];

%% === 3. Gradient Descent with Escape ===
maxSteps = 5000;
tol = 2;
escape_counter = 0;

path_i = start_i;
path_j = start_j;

path = [path_i, path_j];

for step = 1:maxSteps
    
    % Check goal
    if sqrt((path_i - goal_i)^2 + (path_j - goal_j)^2) < tol
        disp('Reached goal!');
        break;
    end
    
    % Find neighbors
    neighbors = [];
    for di = -1:1
        for dj = -1:1
            if di == 0 && dj == 0
                continue;
            end
            ni = path_i + di;
            nj = path_j + dj;
            if ni < 1 || ni > nRows || nj < 1 || nj > nCols
                continue;
            end
            neighbors = [neighbors; ni, nj, U(ni, nj)];
        end
    end
    
    % Pick min potential neighbor
    [~, idx] = min(neighbors(:,3));
    next_i = neighbors(idx,1);
    next_j = neighbors(idx,2);
    
    % If stuck, try random escape
    if next_i == path_i && next_j == path_j
        escape_counter = escape_counter + 1;
        if escape_counter > 10
            % Random move in 8 directions
            ri = randi([-1,1]);
            rj = randi([-1,1]);
            if ~(ri == 0 && rj == 0)
                ni = path_i + ri;
                nj = path_j + rj;
                if ni >= 1 && ni <= nRows && nj >= 1 && nj <= nCols
                    next_i = ni;
                    next_j = nj;
                end
            end
            disp('Escaping local minimum...');
            escape_counter = 0;
        else
            disp('Stuck, trying again...');
        end
    end
    
    % Update
    path_i = next_i;
    path_j = next_j;
    path = [path; path_i, path_j];
end


%% === 4. Plot Path on Potential Field ===
figure;
imagesc(U);
set(gca,'YDir','normal');
colormap(jet); colorbar;
title('Gradient Descent Path on Potential Field');
xlabel('Column j');
ylabel('Row i');
hold on;

plot(path(:,2), path(:,1), 'w-', 'LineWidth', 2);
plot(start_j, start_i, 'go', 'MarkerSize',10, 'LineWidth',2);
plot(goal_j, goal_i, 'rx', 'MarkerSize',12, 'LineWidth',3);
legend('Path','Start','Goal');
