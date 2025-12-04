%% ===================== POTENTIAL FIELD PATH PLANNING ======================
clear; clc; close all;

%% === 1. Load Data =========================================================
risk = readmatrix('riskMap.csv');
lat  = readmatrix('lat.csv');    %#ok<NASGU>
lon  = readmatrix('lon.csv');    %#ok<NASGU>

[nRows, nCols] = size(risk);

%% === 2. User Selects START and GOAL ======================================
figure;
imagesc(risk);
axis equal tight;
set(gca,'YDir','normal');
colormap(jet); colorbar;
title('Click START then GOAL');
xlabel('Column j'); ylabel('Row i');
hold on;

drawnow; shg;

disp('Click START...');
[start_j, start_i] = ginput(1);
start_i = round(start_i); 
start_j = round(start_j);
plot(start_j, start_i, 'go','MarkerSize',12,'LineWidth',2);

disp('Click GOAL...');
[goal_j, goal_i] = ginput(1);
goal_i = round(goal_i); 
goal_j = round(goal_j);
plot(goal_j, goal_i, 'rx','MarkerSize',12,'LineWidth',2);

fprintf('START = (%d,%d), GOAL = (%d,%d)\n', start_i,start_j, goal_i,goal_j);


%% === 3. PROBABILITY-BASED TLS Threshold =================================
% Risk map represents fatality probability per cell.
% TLS sets the maximum acceptable fatality probability.
TLS_thr = 1e-5;   % recommended safety level
fprintf('Using probability-based TLS = %.1e\n', TLS_thr);


%% === 4. Run Potential Field ==============================================
[path_pf, time_pf, metrics_pf] = runPotentialField( ...
    risk, start_i, start_j, goal_i, goal_j, TLS_thr);


%% === 5. Plot PF Path Over Full Map =======================================
figure;
imagesc(risk);
set(gca,'YDir','normal');
axis equal;
colormap(jet); 
colorbar;
colorbar;
caxis([0 1e-5]);
hold on;
plot(path_pf(:,2), path_pf(:,1), 'c-', 'LineWidth',2);
plot(start_j,start_i,'go','MarkerSize',12,'LineWidth',2);
plot(goal_j,goal_i,'rx','MarkerSize',12,'LineWidth',2);

xlim([1 nCols]); 
ylim([1 nRows]);
title(sprintf('Potential Field Path (TLS = %.1e)', TLS_thr));
xlabel('Column j'); ylabel('Row i');
legend('Potential Field Path','Start','Goal');

fprintf('\n===== Potential Field Metrics =====\n');
fprintf('Time (s):          %.4f\n', time_pf);
fprintf('Path Length:       %.2f\n', metrics_pf.length);
fprintf('Max Risk:          %.3g\n', metrics_pf.maxRisk);
fprintf('Mean Risk:         %.3g\n', metrics_pf.meanRisk);
fprintf('Total Fatality P:  %.3g\n', metrics_pf.totalFatalityRisk);
fprintf('TLS Feasible:      %s\n', yesno(metrics_pf.TLS_feasible));


%% ======================= Helper Functions ================================

function out = yesno(x)
    if x, out = 'Yes'; else, out = 'No'; end
end

function [path, time_pf, metrics] = runPotentialField(risk, start_i, start_j, goal_i, goal_j, TLS_thr)

    [nRows,nCols] = size(risk);
    [jGrid,iGrid] = meshgrid(1:nCols, 1:nRows);

% Weights controlling distance vs. risk
w_dist = 1;      % Higher → prefer shorter path
w_risk = 50;     % Higher → avoid risk more strongly

% Distance-to-goal cost
U_att = sqrt((iGrid - goal_i).^2 + (jGrid - goal_j).^2);

% Risk-based repulsive cost (inverse risk → high barrier for danger)
U_rep = w_risk ./ (risk + 1e-12);

% Cap extreme values to avoid numerical walls
U_rep(U_rep > 1e6) = 1e6;

% Combined cost field
U = w_dist * U_att + U_rep;

% Normalize cost field
U = (U - min(U(:))) / (max(U(:)) - min(U(:)) + eps);


    %% === Gradient Descent Search ========================================
    path_i = start_i;
    path_j = start_j;

    prev_i = start_i;
    prev_j = start_j;

    path = [path_i, path_j];
    maxSteps = 8000;
    tol = 2;

    escape_counter = 0;

    tic;
    for step = 1:maxSteps

        % If close to goal, snap to it
        distGoal = hypot(path_i - goal_i, path_j - goal_j);
        if distGoal < 30
            di = sign(goal_i - path_i);
            dj = sign(goal_j - path_j);
            path_i = max(1, min(nRows, path_i + di));
            path_j = max(1, min(nCols, path_j + dj));
            path = [path; path_i, path_j]; %#ok<AGROW>
            if distGoal < tol
                break;
            end
            continue;
        end

        % Check neighbors and move to lowest potential
        bestU = inf;
        next_i = path_i;
        next_j = path_j;

        for di = -1:1
            for dj = -1:1
                if di==0 && dj==0, continue; end
                ni = path_i + di;
                nj = path_j + dj;
                if ni>=1 && ni<=nRows && nj>=1 && nj<=nCols
                    if U(ni,nj) < bestU
                        bestU = U(ni,nj);
                        next_i = ni;
                        next_j = nj;
                    end
                end
            end
        end

        % Escape local minima
        if next_i == path_i && next_j == path_j
            escape_counter = escape_counter + 1;
            if escape_counter > 8
                ri = randi([-1 1]);
                rj = randi([-1 1]);
                ni = path_i + ri;
                nj = path_j + rj;
                if ni>=1 && ni<=nRows && nj>=1 && nj<=nCols
                    next_i = ni; next_j = nj;
                end
                escape_counter = 0;
            end
        end

        % Momentum term
        path_i = round(max(1, min(nRows, next_i + 0.2*(next_i - prev_i))));
        path_j = round(max(1, min(nCols, next_j + 0.2*(next_j - prev_j))));
        prev_i = next_i; prev_j = next_j;

        path = [path; path_i, path_j]; %#ok<AGROW>
    end
    time_pf = toc;

    %% === Evaluate Path ===================================================
    metrics = evaluate_path(path, risk, TLS_thr);
end

function metrics = evaluate_path_MDPI(path, risk, TLS_total)

    pi = path(:,1);
    pj = path(:,2);

    segLen = hypot(diff(pi), diff(pj));
    metrics.length = sum(segLen);

    risk_vals = risk(sub2ind(size(risk), pi, pj));

    metrics.maxRisk  = max(risk_vals);
    metrics.meanRisk = mean(risk_vals);

    % accumulated fatality probability
    metrics.totalFatalityRisk = 1 - prod(1 - risk_vals);

    % MDPI-style mission-level TLS constraint
    metrics.TLS_feasible = metrics.totalFatalityRisk <= TLS_total;

end


