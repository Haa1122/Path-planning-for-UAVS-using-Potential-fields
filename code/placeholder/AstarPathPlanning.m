function AstarPathPlanning
    clear; clc; close all;

    %% 1. Load data
    risk = readmatrix('riskMap.csv');
    [nRows, nCols] = size(risk);

    %% 2. Show map & get clicks
    figure; 
    imagesc(risk); 
    set(gca,'YDir','normal'); 
    colormap(jet); colorbar;
    title('A*: Click START then GOAL');
    xlabel('Column j'); ylabel('Row i');
    hold on;

    disp('Click START point...');
    [sj, si] = ginput(1);
    si = round(si); sj = round(sj);
    plot(sj, si, 'go', 'MarkerSize',12,'LineWidth',2);

    disp('Click GOAL point...');
    [gj, gi] = ginput(1);
    gi = round(gi); gj = round(gj);
    plot(gj, gi, 'rx', 'MarkerSize',12,'LineWidth',2);

    fprintf("Start = (%d,%d), Goal = (%d,%d)\n", si, sj, gi, gj);

    %% 3. A* INITIALIZATION
    open = [];   % [i, j, g, f]
    closed = zeros(nRows, nCols);

    gScore = inf(nRows, nCols);
    gScore(si, sj) = 0;

    fScore = inf(nRows, nCols);
    fScore(si, sj) = heuristic(si, sj, gi, gj);

    open = [si, sj, 0, fScore(si,sj)];

    parent = zeros(nRows, nCols, 2);

    %% 4. A* SEARCH
    moves = [-1 0; 1 0; 0 -1; 0 1; -1 -1; -1 1; 1 -1; 1 1];

    while ~isempty(open)
        [~, idx] = min(open(:,4));
        current = open(idx,:);
        open(idx,:) = [];

        i = current(1);
        j = current(2);

        if i == gi && j == gj
            disp("Goal reached!");
            break;
        end

        closed(i,j) = 1;

        for k = 1:size(moves,1)
            ni = i + moves(k,1);
            nj = j + moves(k,2);

            if ni < 1 || ni > nRows || nj < 1 || nj > nCols
                continue;
            end

            if closed(ni,nj)
                continue;
            end

            % Movement cost with risk
            moveCost = 1 + risk(ni,nj);

            tentative_g = gScore(i,j) + moveCost;

            if tentative_g < gScore(ni,nj)
                parent(ni,nj,:) = [i, j];
                gScore(ni,nj) = tentative_g;
                fScore(ni,nj) = tentative_g + heuristic(ni,nj,gi,gj);

                open = [open; ni, nj, tentative_g, fScore(ni,nj)];
            end
        end
    end

    %% 5. Reconstruct path
    path = [];
    ci = gi; cj = gj;
    while ~(ci == si && cj == sj)
        path = [path; ci, cj];
        p = squeeze(parent(ci,cj,:));
        if p(1) == 0
            disp("No valid path reconstructed!");
            return;
        end
        ci = p(1); cj = p(2);
    end
    path = [path; si, sj];

    %% 6. Plot final path
    figure;
    imagesc(risk);
    set(gca,'YDir','normal');
    colormap(jet); colorbar;
    title('A* Final Path');
    hold on;

    plot(path(:,2), path(:,1), 'w-', 'LineWidth', 2);
    plot(sj, si, 'go', 'MarkerSize',12,'LineWidth',2);
    plot(gj, gi, 'rx', 'MarkerSize',12,'LineWidth',2);
end

%% Heuristic function
function h = heuristic(i, j, gi, gj)
    h = sqrt((i - gi)^2 + (j - gj)^2);
end
