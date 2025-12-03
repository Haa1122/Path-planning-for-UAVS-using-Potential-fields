%% Combined Comparison: Potential Field vs A*
clear; clc; close all;

%% === 1. Load data ===
risk = readmatrix('riskMap.csv');
lat  = readmatrix('lat.csv');   %#ok<NASGU>
lon  = readmatrix('lon.csv');   %#ok<NASGU>

[nRows, nCols] = size(risk);

%% === 2. Click START and GOAL ===
figure;
imagesc(risk);
set(gca,'YDir','normal');
colormap(jet); colorbar;
title('Click START then GOAL');
xlabel('Column j'); ylabel('Row i');
hold on;

drawnow;
shg;

disp('Click START...');
[start_j, start_i] = ginput(1);
start_i = round(start_i); start_j = round(start_j);
plot(start_j, start_i, 'go', 'MarkerSize',12,'LineWidth',2);

disp('Click GOAL...');
[goal_j, goal_i] = ginput(1);
goal_i = round(goal_i); goal_j = round(goal_j);
plot(goal_j, goal_i, 'rx','MarkerSize',12,'LineWidth',2);

fprintf('Start = (%d,%d), Goal = (%d,%d)\n', start_i,start_j, goal_i,goal_j);

%% === 3. TLS Threshold (SAME for both planners) ===
% Relative TLS based on risk-map statistics.
% You can tune this percentile: 85 = stricter, 90 = looser.
TLS_percentile = 85;
TLS_thr = prctile(risk(:), TLS_percentile);

fprintf('TLS percentile = P%d\n', TLS_percentile);
fprintf('TLS threshold  = %.6g\n', TLS_thr);

%% === 4. Run Potential Field ===================================
[path_pf, time_pf, metrics_pf] = runPotentialField( ...
    risk, start_i, start_j, goal_i, goal_j, TLS_thr);

%% === 5. Run A* (TLS-ENFORCED) ================================
% A* parameters (distance vs risk weighting)
w_dist = 1.0;
w_risk = 1.0;

% Enforce TLS as a HARD constraint in A*:
riskBlockThresh = TLS_thr;

tic;
[path_astar, ~] = astar_risk_grid( ...
    risk, [start_i,start_j], [goal_i,goal_j], ...
    w_dist, w_risk, riskBlockThresh);
time_astar = toc;

if isempty(path_astar)
    error('A* did not find a path (under TLS constraint).');
end

metrics_astar = evaluate_path(path_astar, risk, TLS_thr);

%% === 6. Print Comparison Table ===============================
fprintf('\n----------------------------------------------\n');
fprintf('        COMPARISON: POTENTIAL FIELD vs A*\n');
fprintf('          (TLS threshold = %.6g)\n', TLS_thr);
fprintf('----------------------------------------------\n');
fprintf('%-22s %-12s %-12s\n', 'Metric', 'PF', 'A*');
fprintf('----------------------------------------------\n');
fprintf('%-22s %-12.4f %-12.4f\n', 'Time (s)', time_pf, time_astar);
fprintf('%-22s %-12.2f %-12.2f\n', 'Path Length', metrics_pf.length, metrics_astar.length);
fprintf('%-22s %-12.4g %-12.4g\n', 'Accumulated Risk', metrics_pf.accRisk, metrics_astar.accRisk);
fprintf('%-22s %-12.4g %-12.4g\n', 'Max Risk', metrics_pf.maxRisk, metrics_astar.maxRisk);
fprintf('%-22s %-12.4g %-12.4g\n', 'Mean Risk', metrics_pf.meanRisk, metrics_astar.meanRisk);
fprintf('%-22s %-12s %-12s\n', 'TLS Feasible?', yesno(metrics_pf.TLS_feasible), yesno(metrics_astar.TLS_feasible));
fprintf('----------------------------------------------\n');
%% === 7. Plot both paths ======================================
figure;
imagesc(risk);
set(gca,'YDir','normal'); colormap(jet); colorbar;
title('Potential Field (cyan) vs A* (white)');
xlabel('Column j'); ylabel('Row i');
hold on;

plot(path_pf(:,2),    path_pf(:,1),    'c-', 'LineWidth',2);
plot(path_astar(:,2), path_astar(:,1), 'w-', 'LineWidth',2);
plot(start_j,start_i,'go','MarkerSize',10,'LineWidth',2);
plot(goal_j,goal_i,'rx','MarkerSize',10,'LineWidth',2);

legend('Potential Field','A*','Start','Goal');


%% ===== Helper Functions =====================================

function out = yesno(x)
    if x
        out = 'Yes';
    else
        out = 'No';
    end
end

function [path, time_pf, metrics] = runPotentialField(risk, start_i, start_j, goal_i, goal_j, TLS_thr)

    [nRows, nCols] = size(risk);
    [jGrid, iGrid] = meshgrid(1:nCols, 1:nRows);

    % === Attractive potential (toward goal) ===
    k_att = 5e-4;
    U_att = k_att * ((iGrid - goal_i).^2 + (jGrid - goal_j).^2);

    % === Repulsive potential (from high-risk areas) ===
    % IMPORTANT: high risk -> high potential (avoid), low risk -> low potential.
    risk_norm = risk / (max(risk(:)) + eps);   % scale to [0,1]
    k_rep = 0.3;
    U_rep = k_rep * risk_norm;

    % Optional: make areas above TLS_thr even more "repulsive"
    U_rep(risk > TLS_thr) = U_rep(risk > TLS_thr) + 0.5;

    % Total potential
    U = U_att + U_rep;
    U = (U - min(U(:))) / (max(U(:)) - min(U(:)) + eps);

    % Gradient-descent style search
    maxSteps = 8000;
    tol = 2;

    path_i = start_i;
    path_j = start_j;

    prev_i = start_i;
    prev_j = start_j;

    path = [path_i, path_j];
    escape_counter = 0;

    tic;
    for step = 1:maxSteps

        % If close to goal, "snap" toward it
        distGoal = hypot(path_i - goal_i, path_j - goal_j);
        if distGoal < 30
            di = sign(goal_i - path_i);
            dj = sign(goal_j - path_j);
            path_i = max(1,min(nRows, path_i + di));
            path_j = max(1,min(nCols, path_j + dj));
            path = [path; path_i, path_j]; %#ok<AGROW>
            if distGoal < tol
                break;
            end
            continue;
        end

        % Check 8-connected neighbors for lowest potential
        neighbors = [];
        for di = -1:1
            for dj = -1:1
                if di==0 && dj==0
                    continue;
                end
                ni = path_i + di;
                nj = path_j + dj;
                if ni>=1 && ni<=nRows && nj>=1 && nj<=nCols
                    neighbors = [neighbors; ni, nj, U(ni,nj)]; %#ok<AGROW>
                end
            end
        end

        [~, idx] = min(neighbors(:,3));
        next_i = neighbors(idx,1);
        next_j = neighbors(idx,2);

        % Escape local minima
        if next_i==path_i && next_j==path_j
            escape_counter = escape_counter + 1;
            if escape_counter > 8
                ri = randi([-1,1]); rj = randi([-1,1]);
                ni = path_i + ri; nj = path_j + rj;
                if ni>=1 && ni<=nRows && nj>=1 && nj<=nCols
                    next_i = ni; next_j = nj;
                end
                escape_counter = 0;
            end
        end

        % Add a bit of momentum
        path_i = round(max(1,min(nRows, next_i + 0.2*(next_i - prev_i))));
        path_j = round(max(1,min(nCols, next_j + 0.2*(next_j - prev_j))));
        prev_i = next_i; prev_j = next_j;

        path = [path; path_i, path_j]; %#ok<AGROW>
    end
    time_pf = toc;

    metrics = evaluate_path(path, risk, TLS_thr);
end

function metrics = evaluate_path(path, risk, TLS_thr)

    pi = path(:,1);
    pj = path(:,2);

    segLen = hypot(diff(pi), diff(pj));
    metrics.length = sum(segLen);

    risk_vals = risk(sub2ind(size(risk), pi, pj));
    metrics.maxRisk  = max(risk_vals);
    metrics.meanRisk = mean(risk_vals);

    if numel(pi) > 1
        risk_mid = 0.5 * (risk_vals(1:end-1) + risk_vals(2:end));
        metrics.accRisk = sum(segLen .* risk_mid);
    else
        metrics.accRisk = 0;
    end

    metrics.TLS_feasible = all(risk_vals <= TLS_thr);
end

function [path, gScoreMap] = astar_risk_grid(risk, startIJ, goalIJ, w_dist, w_risk, riskBlockThresh)
    % A* on a 2D grid with distance + risk cost, and optional hard blocking.

    [nRows,nCols] = size(risk);
    nCells = nRows * nCols;

    start_i = max(1, min(nRows, startIJ(1)));
    start_j = max(1, min(nCols, startIJ(2)));
    goal_i  = max(1, min(nRows, goalIJ(1)));
    goal_j  = max(1, min(nCols, goalIJ(2)));

    startIdx = sub2ind([nRows,nCols], start_i, start_j);
    goalIdx  = sub2ind([nRows,nCols], goal_i, goal_j);

    if ~isempty(riskBlockThresh)
        blocked = (risk >= riskBlockThresh);
    else
        blocked = false(size(risk));
    end

    gScore = inf(nCells,1);
    fScore = inf(nCells,1);
    cameFrom = zeros(nCells,1,'uint32');
    openSet = false(nCells,1);
    closedSet = false(nCells,1);

    gScore(startIdx) = 0;
    fScore(startIdx) = heuristic_cost([start_i,start_j], [goal_i,goal_j], w_dist);
    openSet(startIdx) = true;

    [dI, dJ] = meshgrid(-1:1, -1:1);
    dI = dI(:); dJ = dJ(:);
    mask = (dI==0 & dJ==0);
    dI(mask) = []; dJ(mask) = [];

    while any(openSet)
        openIdx = find(openSet);
        [~, k] = min(fScore(openIdx));
        currentIdx = openIdx(k);

        if currentIdx == goalIdx
            path = reconstruct_path(cameFrom, currentIdx, [nRows,nCols]);
            gScoreMap = reshape(gScore, nRows, nCols);
            return;
        end

        openSet(currentIdx) = false;
        closedSet(currentIdx) = true;

        [ci, cj] = ind2sub([nRows,nCols], currentIdx);

        for n = 1:length(dI)
            ni = ci + dI(n);
            nj = cj + dJ(n);

            if ni<1 || ni>nRows || nj<1 || nj>nCols
                continue;
            end

            neighIdx = sub2ind([nRows,nCols], ni, nj);

            if closedSet(neighIdx)
                continue;
            end

            if blocked(ni,nj)
                continue;
            end

            stepDist = hypot(double(dI(n)), double(dJ(n)));
            risk_avg = 0.5 * (risk(ci,cj) + risk(ni,nj));
            stepCost = w_dist*stepDist + w_risk*stepDist*risk_avg;

            tentative = gScore(currentIdx) + stepCost;

            if tentative < gScore(neighIdx)
                cameFrom(neighIdx) = currentIdx;
                gScore(neighIdx) = tentative;
                fScore(neighIdx) = tentative + ...
                    heuristic_cost([ni,nj], [goal_i,goal_j], w_dist);
                openSet(neighIdx) = true;
            end
        end
    end

    % If we get here, no path found under TLS constraint
    path = [];
    gScoreMap = reshape(gScore, nRows, nCols);
end

function h = heuristic_cost(a, b, w_dist)
    h = w_dist * hypot(a(1)-b(1), a(2)-b(2));
end

function path = reconstruct_path(cameFrom, currentIdx, gridSize)
    nRows = gridSize(1); %#ok<NASGU>
    rev = zeros(0,2);
    while currentIdx ~= 0
        [ci, cj] = ind2sub(gridSize, currentIdx);
        rev = [rev; ci, cj]; %#ok<AGROW>
        currentIdx = cameFrom(currentIdx);
    end
    path = flipud(rev);
end
