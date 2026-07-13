
%% Stress relaxation viscoelastic property extraction for ligament/tendon testing
% Requires the Optimization Toolbox (lsqcurvefit).
% Fill in the USER INPUTS block below, then run this file.

%% ============================== USER INPUTS ===============================

DATA_FILE = "/Users/najmagure/Downloads/1st Run/DAQ- Axial Force, … - (Timed).txt";

CSA_MM2 = 5.0;          % cross-sectional area of the tissue, mm^2 (from your other script)
L0_MM = 30.0;           % gauge / ligament length, mm (from your other script)

N_PRONY_TERMS = 2;      % number of exponential decay terms in the fit (1-3 typical)

AUTO_DETECT = true;         % try to find the relaxation window automatically first
CONFIRM_AUTO_DETECT = true; % show the detected window and ask for confirmation
                            % before fitting (type 'manual' to override with sliders)

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

    choice = input(sprintf(['Auto-detected relaxation window shown above.\n' ...
        'Press Enter to accept it, or type ''manual'' to adjust with sliders: ']), 's');
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

[popt, r2, fitted] = fitPronySeries(tRel, stress, N_PRONY_TERMS);
fprintf('Fit R^2 = %.5f\n', r2);
reportProperties(popt, strainLevel, N_PRONY_TERMS);

fig = figure('Name', 'Stress relaxation fit', 'Position', [100 100 900 700]);
subplot(3,1,[1 2]);
plot(tRel, stress, '.', 'MarkerSize', 2, 'Color', [0.27 0.51 0.71]);
hold on;
plot(tRel, fitted, '-', 'LineWidth', 2, 'Color', [0.86 0.08 0.24]);
ylabel('Stress (MPa)');
legend('data', sprintf('Prony fit (n=%d), R^2=%.4f', N_PRONY_TERMS, r2), 'Location', 'best');
title('Stress relaxation fit');

subplot(3,1,3);
plot(tRel, stress - fitted, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.7);
yline(0, 'k-', 'LineWidth', 0.5);
xlabel('Time since hold start (s)'); ylabel('Residual (MPa)');

figPath = fullfile(OUTPUT_DIR, "relaxation_fit.png");
saveas(fig, figPath);
fprintf('Saved fit plot to %s\n', figPath);

paramNames = {'sigma_inf'};
for i = 1:N_PRONY_TERMS
    paramNames{end+1} = sprintf('A%d_MPa', i); %#ok<SAGROW>
    paramNames{end+1} = sprintf('tau%d_s', i); %#ok<SAGROW>
end
resultsTable = table(paramNames(:), popt(:), 'VariableNames', {'parameter', 'value'});
resultsPath = fullfile(OUTPUT_DIR, "fit_parameters.csv");
writetable(resultsTable, resultsPath);
fprintf('Saved fit parameters to %s\n', resultsPath);

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
    t = T.Time_s;
    f = T.Force_N;

    fig = figure('Name', 'Select relaxation window', 'Position', [100 100 950 600]);
    ax = axes('Parent', fig, 'Position', [0.1 0.35 0.85 0.55]);
    plot(ax, t, f, 'Color', [0.27 0.51 0.71], 'LineWidth', 0.7);
    xlabel(ax, 'Time (s)'); ylabel(ax, 'Force (N)');
    title(ax, 'Drag Start/End sliders to bracket the stress-relaxation hold, then click Done');
    yl = ylim(ax);

    hold(ax, 'on');
    hPatch = patch(ax, [t(1) t(end) t(end) t(1)], [yl(1) yl(1) yl(2) yl(2)], ...
        [1 0.65 0], 'FaceAlpha', 0.25, 'EdgeColor', 'none');

    data = struct('tStart', t(1), 'tEnd', t(end), 'done', false);
    guidata(fig, data);

    uicontrol('Parent', fig, 'Style', 'text', 'String', 'Start (s)', ...
        'Units', 'normalized', 'Position', [0.05 0.20 0.1 0.04]);
    sStart = uicontrol('Parent', fig, 'Style', 'slider', 'Min', t(1), 'Max', t(end), ...
        'Value', t(1), 'Units', 'normalized', 'Position', [0.16 0.20 0.7 0.04]);

    uicontrol('Parent', fig, 'Style', 'text', 'String', 'End (s)', ...
        'Units', 'normalized', 'Position', [0.05 0.13 0.1 0.04]);
    sEnd = uicontrol('Parent', fig, 'Style', 'slider', 'Min', t(1), 'Max', t(end), ...
        'Value', t(end), 'Units', 'normalized', 'Position', [0.16 0.13 0.7 0.04]);

    uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'Done', ...
        'Units', 'normalized', 'Position', [0.80 0.03 0.12 0.06], ...
        'Callback', @onDone);

    sStart.Callback = @(src, evt) onSliderChange();
    sEnd.Callback = @(src, evt) onSliderChange();

    function onSliderChange()
        lo = min(sStart.Value, sEnd.Value);
        hi = max(sStart.Value, sEnd.Value);
        set(hPatch, 'XData', [lo hi hi lo]);
        drawnow;
    end

    function onDone(~, ~)
        d = guidata(fig);
        d.tStart = min(sStart.Value, sEnd.Value);
        d.tEnd = max(sStart.Value, sEnd.Value);
        d.done = true;
        guidata(fig, d);
        uiresume(fig);
    end

    uiwait(fig);
    d = guidata(fig);
    close(fig);
    if ~d.done
        error('Window closed without clicking Done.');
    end
    window = [d.tStart, d.tEnd];
end

function s = pronyModel(params, t, nTerms)
    sigmaInf = params(1);
    s = sigmaInf * ones(size(t));
    for i = 1:nTerms
        A = params(2*i);
        tau = params(2*i+1);
        s = s + A * exp(-t / tau);
    end
end

function [popt, r2, fitted] = fitPronySeries(tRel, stress, nTerms)
% sigma(t) = sigma_inf + sum(A_i * exp(-t/tau_i)); bounds keep params >= 0.
    sigma0 = max(stress(1), 1e-6);
    tailN = max(5, floor(numel(stress) / 50));
    sigmaInf0 = max(mean(stress(end-tailN+1:end)), 1e-6);
    amp0 = max(sigma0 - sigmaInf0, 1e-6);
    duration = tRel(end);

    tauGuesses = logspace(log10(max(duration/200, 1e-3)), log10(duration/2), nTerms);

    x0 = zeros(1 + 2*nTerms, 1);
    x0(1) = sigmaInf0;
    lb = zeros(1 + 2*nTerms, 1);
    ub = inf(1 + 2*nTerms, 1);
    ub(1) = sigma0;
    for i = 1:nTerms
        x0(2*i) = amp0 / nTerms;
        x0(2*i+1) = tauGuesses(i);
        lb(2*i+1) = 1e-3;
        ub(2*i) = sigma0 * 2;
        ub(2*i+1) = duration * 50;
    end

    model = @(x, t) pronyModel(x, t, nTerms);
    opts = optimoptions('lsqcurvefit', 'Display', 'off', ...
        'MaxFunctionEvaluations', 2e4, 'MaxIterations', 2e4);
    popt = lsqcurvefit(model, x0, tRel, stress, lb, ub, opts);

    fitted = model(popt, tRel);
    ssRes = sum((stress - fitted).^2);
    ssTot = sum((stress - mean(stress)).^2);
    r2 = 1 - ssRes / ssTot;
end

function reportProperties(popt, strainLevel, nTerms)
    sigmaInf = popt(1);
    sigma0 = sigmaInf + sum(popt(2:2:2*nTerms));
    e0 = sigma0 / strainLevel;
    eInf = sigmaInf / strainLevel;
    if sigma0 ~= 0
        pctRelax = (sigma0 - sigmaInf) / sigma0 * 100;
    else
        pctRelax = NaN;
    end

    fprintf('\n--- Fitted viscoelastic properties ---\n');
    fprintf('Equilibrium stress   sigma_inf = %.4f MPa\n', sigmaInf);
    fprintf('Peak (t=0) stress    sigma_0   = %.4f MPa\n', sigma0);
    fprintf('Instantaneous modulus E0       = %.2f MPa\n', e0);
    fprintf('Equilibrium modulus   E_inf    = %.2f MPa\n', eInf);
    fprintf('Percent relaxation             = %.1f %%\n', pctRelax);
    for i = 1:nTerms
        A = popt(2*i); tau = popt(2*i+1);
        fprintf('Term %d: A%d = %.4f MPa   tau%d = %.2f s\n', i, i, A, i, tau);
    end
end
