%% Click-Based Path Planning Using Potential Fields (Final Version)
clear; clc; close all;

%% === 1. Load data ===
risk = readmatrix('riskMap.csv');
lat  = readmatrix('lat.csv');
lon  = readmatrix('lon.csv');

[nRows, nCols] = size(risk);

%% === 2. Show risk map and GET user clicks ===
figure;
imagesc(risk);     
set(gca,'YDir','normal');
colormap(jet); colorbar;
title('Click START then GOAL');
xlabel('Column j'); ylabel('Row i');
hold on;

disp('Click START point...');
[start_j, start_i] = ginput(1);
start_i = round(start_i); start_j = round(start_j);
plot(start_j, start_i, 'go', 'MarkerSize',12,'LineWidth',2);

disp('Click GOAL point...');
[goal_j, goal_i] = ginput(1);
goal_i = round(goal_i); goal_j = round(goal_j);
plot(goal_j, goal_i, 'rx','MarkerSize',12,'LineWidth',2);

fprintf('Start = (%d,%d), Goal = (%d,%d)\n', start_i,start_j, goal_i,goal_j);

%% === 3. Meshgrid ===
[jGrid, iGrid] = meshgrid(1:nCols, 1:nRows);

%% === 4. Attractive potential ===
k_att = 5e-4;   % stronger attraction
U_att = k_att * ((iGrid - goal_i).^2 + (jGrid - goal_j).^2);

%% === 5. Repulsive potential using softened log-risk ===
risk_clipped = max(risk, 1e-12);
risk_soft = risk_clipped + 5e-3;     % smoother near goal
risk_log = -log(risk_soft);          % repulsion
risk_log = risk_log / max(risk_log(:));

k_rep = 0.3;                         % balanced repulsion
U_rep = k_rep * risk_log;

%% === 6. Total potential (normalize) ===
U = U_att + U_rep;
U = (U - min(U(:))) / (max(U(:)) - min(U(:)));

%% DEBUG (optional)
figure; imagesc(U);
set(gca,'YDir','normal'); colorbar;
title('DEBUG: Potential Field U');
pause(0.2);

%% === 7. Gradient Descent Path Planning ===
maxSteps = 8000;
tol = 2;

path_i = start_i;
path_j = start_j;

prev_i = start_i;
prev_j = start_j;

path = [path_i, path_j];
escape_counter = 0;

for step = 1:maxSteps

    % === Goal magnet zone (forces completion) ===
    distGoal = sqrt((path_i - goal_i)^2 + (path_j - goal_j)^2);
    if distGoal < 30
        di = sign(goal_i - path_i);
        dj = sign(goal_j - path_j);

        next_i = path_i + di;
        next_j = path_j + dj;

        % clamp
        next_i = max(1, min(nRows, next_i));
        next_j = max(1, min(nCols, next_j));

        path = [path; next_i, next_j];
        path_i = next_i; 
        path_j = next_j;

        if distGoal < tol
            disp("Reached goal!");
            break;
        end

        continue;
    end

    % === Normal gradient descent ===
    neighbors = [];
    for di = -1:1
        for dj = -1:1
            if di == 0 && dj == 0, continue; end
            ni = path_i + di; nj = path_j + dj;
            if ni < 1 || ni > nRows || nj < 1 || nj > nCols
                continue;
            end
            neighbors = [neighbors; ni, nj, U(ni,nj)];
        end
    end

    [~, idx] = min(neighbors(:,3));
    next_i = neighbors(idx,1);
    next_j = neighbors(idx,2);

    % Escape from local trap
    if next_i == path_i && next_j == path_j
        escape_counter = escape_counter + 1;
        if escape_counter > 8
            ri = randi([-1,1]); 
            rj = randi([-1,1]);
            if ~(ri==0 && rj==0)
                ni = path_i + ri; nj = path_j + rj;
                if ni>=1 && ni<=nRows && nj>=1 && nj<=nCols
                    next_i = ni; next_j = nj;
                end
            end
            escape_counter = 0;
        end
    end

    % momentum
    path_i = next_i + 0.2*(next_i - prev_i);
    path_j = next_j + 0.2*(next_j - prev_j);

    prev_i = next_i; prev_j = next_j;

    % clamp
    path_i = round(max(1, min(nRows, path_i)));
    path_j = round(max(1, min(nCols, path_j)));

    path = [path; path_i, path_j];
end

%% === 8. Plot final path ===
figure;
imagesc(risk);
set(gca,'YDir','normal'); colormap(jet); colorbar;
title('Final Potential-Field Path');
xlabel('Column j'); ylabel('Row i');
hold on;

plot(path(:,2), path(:,1), 'w-', 'LineWidth',2);
plot(start_j,start_i,'go','MarkerSize',12,'LineWidth',2);
plot(goal_j,goal_i,'rx','MarkerSize',12,'LineWidth',2);
legend('Path','Start','Goal');

