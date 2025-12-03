%% ===================== STRICT-TLS POTENTIAL FIELD ======================
clear; clc; close all;

%% 1. Load data
risk = readmatrix('riskMap.csv');
lat  = readmatrix('lat.csv'); %#ok<NASGU>
lon  = readmatrix('lon.csv'); %#ok<NASGU>

[nRows,nCols] = size(risk);

%% 2. TLS: strict fatality probability threshold
TLS_thr = 1e-5;   % DO NOT CHANGE if you want strict TLS
fprintf('Using strict probability-based TLS = %.1e\n', TLS_thr);

%% 3. Pick START and GOAL (user clicks)
figure;
imagesc(risk);
set(gca,'YDir','normal'); axis equal tight;
colormap(jet); colorbar;
caxis([0 TLS_thr]);                 % scale 0 → 10^-5 as you requested
title('Click START then GOAL (display clipped at 10^{-5})');
xlabel('Column j'); ylabel('Row i');
hold on; drawnow;

disp('Click START...');
[start_j,start_i] = ginput(1);
start_i = round(start_i); start_j = round(start_j);
plot(start_j,start_i,'wo','MarkerSize',10,'LineWidth',2);

disp('Click GOAL...');
[goal_j,goal_i] = ginput(1);
goal_i = round(goal_i); goal_j = round(goal_j);
plot(goal_j,goal_i,'wx','MarkerSize',10,'LineWidth',2);

fprintf('Raw START = (%d,%d), raw GOAL = (%d,%d)\n', start_i,start_j, goal_i,goal_j);

%% 4. Move START/GOAL to nearest TLS-safe cells if necessary
[start_i,start_j, movedS] = moveToNearestSafe(risk, start_i,start_j, TLS_thr);
[goal_i, goal_j, movedG] = moveToNearestSafe(risk, goal_i, goal_j, TLS_thr);

if movedS
    fprintf('START was unsafe; moved to nearest TLS-safe cell (%d,%d).\n', start_i,start_j);
end
if movedG
    fprintf('GOAL was unsafe; moved to nearest TLS-safe cell (%d,%d).\n', goal_i,goal_j);
end

plot(start_j,start_i,'go','MarkerSize',12,'LineWidth',2); % safe start
plot(goal_j,goal_i,'rx','MarkerSize',12,'LineWidth',2);   % safe goal
legend('Raw Start','Raw Goal','TLS-safe Start','TLS-safe Goal','Location','best');

%% 5. Run TLS-strict potential field
[path_pf, time_pf, metrics_pf] = runPotentialField_TLSstrict( ...
    risk, start_i, start_j, goal_i, goal_j, TLS_thr);

if isempty(path_pf)
    fprintf('\nNo TLS-safe path exists between START and GOAL under TLS = %.1e\n', TLS_thr);
    return;
end

%% 6. Plot final path with same 0→1e-5 scale
figure;
imagesc(risk);
set(gca,'YDir','normal'); axis equal;
colormap(jet); colorbar;
caxis([0 TLS_thr]);        % keep scale at 10^-5
hold on;

xlim([1 nCols]); ylim([1 nRows]);
plot(path_pf(:,2), path_pf(:,1),'c-','LineWidth',2);
plot(start_j,start_i,'go','MarkerSize',12,'LineWidth',2);
plot(goal_j,goal_i,'rx','MarkerSize',12,'LineWidth',2);

title(sprintf('TLS-Strict Potential Field Path (TLS = %.1e)', TLS_thr));
xlabel('Column j'); ylabel('Row i');
legend('Potential Field Path','Start','Goal','Location','best');

%% 7. Print metrics
fprintf('\n===== TLS-STRICT POTENTIAL FIELD METRICS =====\n');
fprintf('Time (s):          %.4f\n', time_pf);
fprintf('Path Length:       %.2f (grid units)\n', metrics_pf.length);
fprintf('Max Risk:          %.3g\n', metrics_pf.maxRisk);
fprintf('Mean Risk:         %.3g\n', metrics_pf.meanRisk);
fprintf('Total Fatality P:  %.3g\n', metrics_pf.totalFatalityRisk);
fprintf('TLS Feasible:      %s\n', yesno(metrics_pf.TLS_feasible));


%% ============== Helper: yes/no ==========================================
function out = yesno(x)
    if x, out = 'Yes'; else, out = 'No'; end
end


%% ============== Helper: move point to nearest TLS-safe cell =============
function [si,sj, moved] = moveToNearestSafe(risk, si,sj, TLS_thr)

[nRows,nCols] = size(risk);
if risk(si,sj) <= TLS_thr
    moved = false;
    return;
end

moved = true;
maxR = max(nRows,nCols);

for r = 1:maxR
    imin = max(1, si-r); imax = min(nRows, si+r);
    jmin = max(1, sj-r); jmax = min(nCols, sj+r);

    sub = risk(imin:imax, jmin:jmax);
    [ii,jj] = find(sub <= TLS_thr);
    if ~isempty(ii)
        gi = ii + imin - 1;
        gj = jj + jmin - 1;
        d  = hypot(gi - si, gj - sj);
        [~,idx] = min(d);
        si = gi(idx);
        sj = gj(idx);
        return;
    end
end

error('No TLS-safe cell found in the entire map (this should not happen).');
end


%% ============== TLS-strict Potential Field ==============================
function [path, time_pf, metrics] = runPotentialField_TLSstrict( ...
    risk, start_i, start_j, goal_i, goal_j, TLS_thr)

[nRows,nCols] = size(risk);
[jGrid,iGrid] = meshgrid(1:nCols, 1:nRows);

% Attractive potential
k_att = 5e-4;
U_att = k_att*((iGrid-goal_i).^2 + (jGrid-goal_j).^2);

% Repulsive potential (weak, so drone can pass near high risk)
risk_norm = risk / (max(risk(:)) + eps);
k_rep = 0.02;
U_rep = k_rep * risk_norm;

% Optional: cliff at TLS to shape gradients (but movement check enforces safety)
U_rep(risk > TLS_thr) = U_rep(risk > TLS_thr) + 5.0;

% Total potential
U = U_att + U_rep;
U = (U - min(U(:))) / (max(U(:))-min(U(:)) + eps);

% PF search using "best safe neighbor" descent
path_i = start_i;
path_j = start_j;
path   = [path_i, path_j];

maxSteps = 15000;
tol = 2;

tic;
for step = 1:maxSteps

    % Check if close enough to goal
    if hypot(path_i - goal_i, path_j - goal_j) < tol
        break;
    end

    bestU   = Inf;
    next_i  = path_i;
    next_j  = path_j;

    % examine 8-connected neighbors
    for di = -1:1
        for dj = -1:1
            if di==0 && dj==0, continue; end
            ni = path_i + di;
            nj = path_j + dj;

            if ni < 1 || ni > nRows || nj < 1 || nj > nCols
                continue;
            end

            % STRICT TLS: never enter unsafe cells
            if risk(ni,nj) > TLS_thr
                continue;
            end

            if U(ni,nj) < bestU
                bestU = U(ni,nj);
                next_i = ni;
                next_j = nj;
            end
        end
    end

    % if we couldn't find any better safe neighbor, stop
    if next_i == path_i && next_j == path_j
        % no safe descent direction → local minimum or dead end
        break;
    end

    path_i = next_i;
    path_j = next_j;
    path   = [path; path_i, path_j]; %#ok<AGROW>
end
time_pf = toc;

% If we did not get close to goal, declare failure
if hypot(path_i - goal_i, path_j - goal_j) >= tol
    path = [];
    metrics = struct('length',0,'maxRisk',NaN,'meanRisk',NaN, ...
                     'totalFatalityRisk',NaN,'TLS_feasible',false);
    return;
end

% Evaluate path
metrics = evaluate_path_TLS(risk, path, TLS_thr);

end


%% ============== Path Evaluation (TLS-safe) ==============================
function metrics = evaluate_path_TLS(risk, path, TLS_thr)

pi = path(:,1);
pj = path(:,2);

segLen = hypot(diff(pi), diff(pj));
metrics.length = sum(segLen);

risk_vals = risk(sub2ind(size(risk), pi, pj));
metrics.maxRisk  = max(risk_vals);
metrics.meanRisk = mean(risk_vals);

metrics.totalFatalityRisk = 1 - prod(1 - risk_vals);

metrics.TLS_feasible = all(risk_vals <= TLS_thr);

end
