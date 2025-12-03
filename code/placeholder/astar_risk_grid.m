function [path, gScoreMap] = astar_risk_grid(risk, startIJ, goalIJ, w_dist, w_risk, riskBlockThresh)
% A* search on a 2D risk grid with 8-connected neighbors.
%
% INPUTS:
%   risk   : [nRows x nCols] ground-risk map (non-negative)
%   startIJ: [1x2] [i j] start indices
%   goalIJ : [1x2] [i j] goal indices
%   w_dist : weight for distance in the cost function
%   w_risk : weight for risk in the cost function
%   riskBlockThresh: cells with risk >= this are treated as blocked
%                    (use [] to disable)
%
% OUTPUTS:
%   path      : [k x 2] matrix of [i j] indices from start to goal
%   gScoreMap : [nRows x nCols] map of g-scores (cost from start)
%
% Cost per step:
%   step_cost = w_dist * step_distance + w_risk * step_distance * risk_avg

    if nargin < 6
        riskBlockThresh = [];
    end

    [nRows, nCols] = size(risk);
    nCells = nRows * nCols;

    % Clamp start/goal inside grid
    start_i = max(1, min(nRows, startIJ(1)));
    start_j = max(1, min(nCols, startIJ(2)));
    goal_i  = max(1, min(nRows, goalIJ(1)));
    goal_j  = max(1, min(nCols, goalIJ(2)));

    startIdx = sub2ind([nRows,nCols], start_i, start_j);
    goalIdx  = sub2ind([nRows,nCols], goal_i,  goal_j);

    % Precompute blocked cells
    if ~isempty(riskBlockThresh)
        blocked = risk >= riskBlockThresh;
    else
        blocked = false(size(risk));
    end

    % A* arrays
    gScore = inf(nCells,1);
    fScore = inf(nCells,1);
    cameFrom = zeros(nCells,1,'uint32');
    openSet = false(nCells,1);
    closedSet = false(nCells,1);

    gScore(startIdx) = 0;
    fScore(startIdx) = heuristic_cost([start_i,start_j],[goal_i,goal_j],w_dist);
    openSet(startIdx) = true;

    % Neighbor offsets (8-connected)
    [dI, dJ] = meshgrid(-1:1, -1:1);
    dI = dI(:); dJ = dJ(:);
    zeroMask = (dI==0 & dJ==0);
    dI(zeroMask) = []; dJ(zeroMask) = [];

    while any(openSet)
        % Get node in openSet with lowest fScore
        openIdx = find(openSet);
        [~,k] = min(fScore(openIdx));
        currentIdx = openIdx(k);

        if currentIdx == goalIdx
            % reconstruct path
            path = reconstruct_path(cameFrom, currentIdx, [nRows,nCols]);
            gScoreMap = reshape(gScore, nRows, nCols);
            return;
        end

        openSet(currentIdx) = false;
        closedSet(currentIdx) = true;

        [ci, cj] = ind2sub([nRows,nCols], currentIdx);

        % Explore neighbors
        for n = 1:length(dI)
            ni = ci + dI(n);
            nj = cj + dJ(n);

            if ni < 1 || ni > nRows || nj < 1 || nj > nCols
                continue;
            end

            neighIdx = sub2ind([nRows,nCols], ni, nj);

            if closedSet(neighIdx)
                continue;
            end

            if blocked(ni,nj)
                continue;  % treated as obstacle
            end

            % Cost from current to neighbor
            stepDist = hypot(double(dI(n)), double(dJ(n)));

            risk_curr = risk(ci,cj);
            risk_neigh = risk(ni,nj);
            risk_avg = 0.5 * (risk_curr + risk_neigh);

            stepCost = w_dist * stepDist + w_risk * stepDist * risk_avg;

            tentative_gScore = gScore(currentIdx) + stepCost;

            if tentative_gScore < gScore(neighIdx)
                cameFrom(neighIdx) = currentIdx;
                gScore(neighIdx) = tentative_gScore;
                fScore(neighIdx) = tentative_gScore + ...
                    heuristic_cost([ni,nj],[goal_i,goal_j],w_dist);

                if ~openSet(neighIdx)
                    openSet(neighIdx) = true;
                end
            end
        end
    end

    % If we reach here, no path found
    path = [];
    gScoreMap = reshape(gScore, nRows, nCols);
end

function h = heuristic_cost(nodeIJ, goalIJ, w_dist)
    % Straight-line distance heuristic (scaled by w_dist)
    di = double(nodeIJ(1) - goalIJ(1));
    dj = double(nodeIJ(2) - goalIJ(2));
    h = w_dist * hypot(di, dj);
end

function path = reconstruct_path(cameFrom, currentIdx, gridSize)
    nRows = gridSize(1);
    % Reconstruct in reverse order
    revPath = zeros(0,2);
    while currentIdx ~= 0
        [ci,cj] = ind2sub([nRows,gridSize(2)], currentIdx);
        revPath = [revPath; ci, cj]; %#ok<AGROW>
        currentIdx = cameFrom(currentIdx);
    end
    path = flipud(revPath);
end