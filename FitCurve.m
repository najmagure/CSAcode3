%% Stress relaxation viscoelastic property extraction
%% ============================== USER INPUTS ===============================

DATA_FILE = "/Users/najmagure/Downloads/1st Run/DAQ- Axial Force, … - (Timed).txt";

CSA_MM2 = 5.0;          % cross-sectional area of the tissue, mm^2 (from your other script)
L0_MM = 30.0;           % gauge / ligament length, mm (from your other script)

AUTO_DETECT = true;         % try to find the relaxation window automatically
CONFIRM_AUTO_DETECT = true; % show the detected window before fitting (press Enter to continue)

OUTPUT_DIR = fullfile(fileparts(mfilename('fullpath')), "output");

if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
end

T = loadDaqFile(DATA_FILE);
T = computeStressStrain(T, CSA_MM2, L0_MM);

window = [];
if AUTO_DETECT
    window = autoDetectRelaxation(T);
end

if ~isempty(window) && CONFIRM_AUTO_DETECT
    fig = figure('Name', 'Auto-detected relaxation window');
    plot(T.Time_s, T.Force_N, 'Color', [0.27 0.51 0.71], 'LineWidth', 0.6);
    hold on;
    yl = ylim;
    patch(window([1 2 2 1]), [yl(1) yl(1) yl(2) yl(2)], ...
        [1 0.65 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    xlabel('Time (s)'); ylabel('Force (N)');
    title(sprintf('Auto-detected window: %.2fs - %.2fs (duration %.1fs)', ...
        window(1), window(2), window(2) - window(1)));
    drawnow;
    choice = input('Press Enter to accept, or type ''manual'' to adjust: ', 's');
    close(fig);
    if strcmpi(strtrim(choice), 'manual')
        window = [];
    end
end

if isempty(window)
    disp('Opening manual selector...');
    window = manualSelectRelaxation(T);
end

tStart = window(1); tEnd = window(2);
mask = T.Time_s >= tStart & T.Time_s <= tEnd;
seg = T(mask, :);
fprintf('\nUsing relaxation window: %.3fs to %.3fs (%d points)\n', tStart, tEnd, height(seg));

tRel = seg.Time_s - seg.Time_s(1);
stress = seg.Stress_MPa;
strainLevel = mean(seg.Strain);

[popt, r2, fitted] = fitQLV(tRel, stress);
fprintf('Fit R^2 = %.5f\n', r2);

resultsT = qlvProperties(popt, strainLevel);
disp(resultsT);

fig = figure('Name', 'Stress relaxation fit', 'Position', [100 100 900 600]);
plot(tRel, stress, '.', 'MarkerSize', 2, 'Color', [0.27 0.51 0.71]);
hold on;
plot(tRel, fitted, '-', 'LineWidth', 2, 'Color', [0.86 0.08 0.24]);
xlabel('Time since hold start (s)'); ylabel('Stress (MPa)');
legend('data', sprintf('QLV fit, R^2=%.4f', r2), 'Location', 'best');
title('Stress relaxation fit');

figPath = fullfile(OUTPUT_DIR, "relaxation_fit.png");
saveas(fig, figPath);
fprintf('Saved fit plot to %s\n', figPath);

resultsPath = fullfile(OUTPUT_DIR, "viscoelastic_properties.csv");
writetable(resultsT, resultsPath);
fprintf('Saved viscoelastic properties to %s\n', resultsPath);

seg.Time_rel_s = tRel;
segPath = fullfile(OUTPUT_DIR, "relaxation_segment.csv");
writetable(seg, segPath);
fprintf('Saved relaxation segment data to %s\n', segPath);


%% ============================== LOCAL FUNCTIONS ============================

function T = loadDaqFile(path)
    fid = fopen(path, 'rb');
    if fid == -1
        error('Could not open data file (check the path/filename exactly): %s', path);
    end
    magic = fread(fid, 2, 'uint8=>char')';
    fclose(fid);

    if strcmp(magic, 'PK')
        raw = readcell(path, 'FileType', 'spreadsheet');
    else
        raw = readcell(path, 'FileType', 'text', 'Delimiter', '\t');
    end

    headerRow = 0;
    for i = 1:min(20, size(raw,1))
        for j = 1:size(raw,2)
            v = raw{i,j};
            if (ischar(v) || isstring(v)) && strcmpi(strtrim(string(v)), "time")
                headerRow = i;
                break;
            end
        end
        if headerRow > 0
            break;
        end
    end
    if headerRow == 0
        error('Could not find a header row with a ''Time'' column in the first 20 rows of %s', path);
    end

    ncols = 0;
    for j = 1:size(raw,2)
        v = raw{headerRow, j};
        if isMissingCell(v)
            break;
        end
        ncols = ncols + 1;
    end

    colNames = strings(1, ncols);
    for j = 1:ncols
        colNames(j) = strtrim(string(raw{headerRow, j}));
    end

    dataRows = raw(headerRow+2:end, 1:ncols);  % +2 skips the units row
    nrows = size(dataRows, 1);
    M = nan(nrows, ncols);
    for r = 1:nrows
        for c = 1:ncols
            v = dataRows{r,c};
            if isnumeric(v) && isscalar(v)
                M(r,c) = v;
            elseif ischar(v) || isstring(v)
                M(r,c) = str2double(v);
            end
        end
    end
    M = M(~any(isnan(M), 2), :);

    if isempty(M)
        safeColNames = colNames;
        safeColNames(ismissing(safeColNames)) = "(blank)";
        error('No numeric data rows parsed from %s (header row %d, columns: %s).', ...
            path, headerRow, strjoin(safeColNames, ', '));
    end

    lowerNames = lower(colNames);
    idxForce = find(contains(lowerNames, "force"), 1);
    idxDisp  = find(contains(lowerNames, "displacement"), 1);
    idxTime  = find(contains(lowerNames, "time"), 1);
    if isempty(idxForce) || isempty(idxDisp) || isempty(idxTime)
        error('Could not identify Force/Displacement/Time columns from header: %s', strjoin(colNames, ', '));
    end

    T = table(M(:,idxTime), M(:,idxForce), M(:,idxDisp), ...
        'VariableNames', {'Time_s', 'Force_N', 'Displacement_mm'});
    T = sortrows(T, 'Time_s');
end

function tf = isMissingCell(v)
% readcell returns MATLAB's dedicated `missing` type for a blank
% spreadsheet cell, not NaN/'' - ismissing() is what catches that.
    tf = isempty(v) || any(ismissing(v)) || ...
        ((ischar(v) || isstring(v)) && strlength(strtrim(string(v))) == 0);
end

function T = computeStressStrain(T, csaMm2, l0Mm)
    T.Stress_MPa = T.Force_N / csaMm2;
    T.Strain = T.Displacement_mm / l0Mm;
end

function window = autoDetectRelaxation(T, minHoldS, minDispFraction, settleTrimS, stillWindowS, stillFrac)
% Longest constant-displacement hold at a non-trivial displacement level.
% Returns [tStart tEnd], or [] if nothing qualifies.
    if nargin < 2, minHoldS = 15.0; end
    if nargin < 3, minDispFraction = 0.3; end
    if nargin < 4, settleTrimS = 0.3; end
    if nargin < 5, stillWindowS = 2.0; end
    if nargin < 6, stillFrac = 0.01; end

    t = T.Time_s;
    d = T.Displacement_mm;

    dt = median(diff(t));
    win = max(3, round(stillWindowS / dt));
    rollStd = movstd(d, win);

    dispRange = max(d) - min(d);
    stillThresh = stillFrac * dispRange;
    isStill = rollStd < stillThresh;

    dispThresh = minDispFraction * max(d);

    n = numel(t);
    bestDuration = -inf;
    bestStart = NaN; bestEnd = NaN;
    i = 1;
    while i <= n
        if isStill(i)
            j = i;
            while j <= n && isStill(j)
                j = j + 1;
            end
            segT0 = t(i); segT1 = t(j-1);
            duration = segT1 - segT0;
            segMeanDisp = mean(d(i:j-1));
            if duration >= minHoldS && segMeanDisp >= dispThresh && duration > bestDuration
                bestDuration = duration;
                bestStart = segT0;
                bestEnd = segT1;
            end
            i = j;
        else
            i = i + 1;
        end
    end

    if isnan(bestStart)
        window = [];
        return;
    end

    bestStart = min(bestStart + settleTrimS, bestEnd);  % trim ramp-settling transient
    window = [bestStart, bestEnd];
end

function window = manualSelectRelaxation(T)
% Draggable drawline ROIs (Image Processing Toolbox), same pattern as
% CalculatingCSA.m's region selector.
    t = T.Time_s;
    f = T.Force_N;

    fig = figure('Name', 'Select relaxation window');
    plot(t, f, 'Color', [0.27 0.51 0.71], 'LineWidth', 0.7);
    hold on;
    xlabel('Time (s)'); ylabel('Force (N)');
    title({'Drag the red lines to bracket the stress-relaxation hold'; 'Click "Done" when finished'});

    tmin = min(t); tmax = max(t);
    t1 = tmin + 0.25 * (tmax - tmin);
    t2 = tmin + 0.75 * (tmax - tmin);

    h1 = drawline('Position', [t1 min(ylim); t1 max(ylim)], 'Color', 'r', 'LineWidth', 2);
    h2 = drawline('Position', [t2 min(ylim); t2 max(ylim)], 'Color', 'r', 'LineWidth', 2);

    uicontrol('Style', 'pushbutton', 'String', 'Done', ...
        'Units', 'normalized', 'Position', [0.45 0.01 0.1 0.06], ...
        'Callback', 'uiresume(gcbf)');
    uiwait(fig);

    window = sort([h1.Position(1,1), h2.Position(1,1)]);
    close(fig);
end

function s = qlvModel(p, t)
% Fung's QLV reduced relaxation function: sigma(t) = sigma0 * G(t).
    sigma0 = p(1); c = p(2); tau1 = p(3); tau2 = p(3) + p(4);
    tSafe = max(t, 1e-6);  % E1(0) diverges; the true limit G(0)=1 is reached as t->0+
    G = (1 + c * (expint(tSafe/tau2) - expint(tSafe/tau1))) / (1 + c*log(tau2/tau1));
    s = sigma0 * G;
end

function [popt, r2, fitted] = fitQLV(tRel, stress)
    sigma0_0 = max(stress(1), 1e-6);
    duration = tRel(end);
    tau1_0 = max(duration/1000, 1e-3);
    tau2_0 = duration * 5;

    x0 = [sigma0_0; 1; tau1_0; tau2_0 - tau1_0];  % [sigma0, c, tau1, tau2-tau1]
    lb = [0; 1e-4; 1e-3; 1e-3];
    ub = [sigma0_0 * 2; 100; duration; duration * 1000];

    opts = optimoptions('lsqcurvefit', 'Display', 'off', ...
        'MaxFunctionEvaluations', 2e4, 'MaxIterations', 2e4);
    popt = lsqcurvefit(@qlvModel, x0, tRel, stress, lb, ub, opts);

    fitted = qlvModel(popt, tRel);
    ssRes = sum((stress - fitted).^2);
    ssTot = sum((stress - mean(stress)).^2);
    r2 = 1 - ssRes / ssTot;
end

function resultsT = qlvProperties(popt, strainLevel)
    c = popt(2); tau1 = popt(3); tau2 = popt(3) + popt(4);
    e0 = popt(1) / strainLevel;
    eInf = e0 / (1 + c * log(tau2 / tau1));
    pctRelax = (1 - eInf / e0) * 100;

    Property = {'Instantaneous modulus (E0)'; 'Equilibrium modulus (E_inf)'; ...
        'Percent relaxation'; 'Damping coefficient (c)'};
    Value = round([e0; eInf; pctRelax; c], 4);
    Units = {'MPa'; 'MPa'; '%'; '-'};
    resultsT = table(Property, Value, Units);
end
