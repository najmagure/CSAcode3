%% DIC_Acumen_Pipeline.m
%
% Automates the manual workflow from "7-1 Data with Plots.xlsx":
%   1. Import raw DIC point-cloud data (exported "name/type/attribute/id" format)
%   2. Import Acumen/MTS DAQ data (force, displacement, time)
%   3. Build "insertion" and "origin" centroids from tracked DIC points
%   4. Compute distance, displacement, and uniaxial strain from the centroids
%   5. Resample DIC and Acumen time series onto a common, user-chosen frequency
%   6. Plot Time vs Strain (DIC + Acumen) and Time vs Force

scriptDir  = fileparts(mfilename("fullpath"));
OUTPUT_DIR = fullfile(scriptDir, "dic_acumen_output");
if ~exist(OUTPUT_DIR, "dir")
    mkdir(OUTPUT_DIR);
end

%% ===================== CONFIG =====================
[dicName, dicPath] = uigetfile({'*.csv;*.xlsx;*.xls', 'DIC data files (*.csv, *.xlsx, *.xls)'; '*.*', 'All files'}, ...
    'Select raw DIC data file', scriptDir);
if isequal(dicName, 0)
    error('No DIC file selected.');
end
cfg.dicFile        = fullfile(dicPath, dicName);   % raw DIC export, as received from the machine
cfg.dicSheet       = 1;                     % sheet name or index (only used if an .xlsx is selected)

% Point IDs (0-indexed, matching the DIC export's "id" column) that make up
% each centroid. Defaults below were reverse-engineered from the manual
% "PCL DIC Data" sheet (Centroid insertion = pts 0,1,2; Centroid origin = pts 6,7,8).
cfg.insertionPointIDs = [0 1 2];
cfg.originPointIDs    = [6 7 8];

% Acumen DAQ file for this run, in the format as originally exported by the
% machine: often a .txt-named file that is actually a zipped .xlsx under the
% hood (no manual conversion needed - the PK signature is auto-detected below).
[acumenName, acumenPath] = uigetfile({'*.txt;*.xlsx;*.xls', 'Acumen DAQ files (*.txt, *.xlsx, *.xls)'; '*.*', 'All files'}, ...
    'Select Acumen DAQ data file', scriptDir);
if isequal(acumenName, 0)
    error('No Acumen file selected.');
end
cfg.acumenFile = fullfile(acumenPath, acumenName);

% The DAQ's own clock rarely starts at the same instant as the DIC video.
% Set acumenTimeOffset (seconds) so that (acumenTime - acumenTimeOffset) lines
% up with DIC time = 0 at the first frame. "auto" zeros to the first DAQ sample.
cfg.acumenTimeOffset = "auto";

cfg.outputXlsx       = fullfile(OUTPUT_DIR, "Resampled_Results.xlsx");

%% ===================== PIPELINE =====================
dic = importDICRaw(cfg.dicFile, cfg.dicSheet);
numFrames = size(dic.X, 2);

acumen = importAcumenDAQ(cfg.acumenFile);
if cfg.acumenTimeOffset == "auto"
    acumenTimeOffset = acumen.Time(1);
else
    acumenTimeOffset = cfg.acumenTimeOffset;
end
acumenTime = acumen.Time - acumenTimeOffset;

% This analysis focuses on the stress-relaxation hold, not the whole test
% (ramp + cyclic loading included). Auto-detect the longest constant-
% displacement hold from the Acumen data, same approach as FitCurve.m's
% autoDetectRelaxation, then show it for a visual sanity check.
relaxWindow = autoDetectRelaxationWindow(acumenTime, acumen.AxialDisplacement);

assumedFs = 30;   % just for turning the detected window into a starting suggestion below
if ~isempty(relaxWindow)
    suggestedRef = max(1, round(relaxWindow(1) * assumedFs) + 1);
    suggestedEnd = min(numFrames, round(relaxWindow(2) * assumedFs) + 1);
    previewFig = figure("Name", "Auto-detected relaxation window", "Color", "white");
    plot(acumenTime, acumen.AxialDisplacement, "Color", [0.27 0.51 0.71], "LineWidth", 0.7);
    hold on
    yl = ylim;
    patch(relaxWindow([1 2 2 1]), [yl(1) yl(1) yl(2) yl(2)], [1 0.65 0], "FaceAlpha", 0.3, "EdgeColor", "none");
    xlabel("Time (s)"); ylabel("Axial Displacement");
    title(sprintf("Auto-detected hold: %.2fs - %.2fs (adjust the frame numbers below if this looks wrong)", ...
        relaxWindow(1), relaxWindow(2)));
    drawnow;
else
    suggestedRef = 1;
    suggestedEnd = numFrames;
    warning("Could not auto-detect a relaxation hold; defaulting to the full test range - set the frame range manually below.");
end

% Reference frame, analysis end frame, and both frequencies are picked
% interactively per test run instead of hardcoded, since these vary test
% to test (matches the picker used in StrainPlot.m). Reference/end frame
% default to the auto-detected relaxation window above, but can be typed
% over manually if the auto-detection picked the wrong section.
answer = inputdlg( ...
    {'Reference Frame (T = 0)', 'Last Frame to Analyze', 'DIC Frequency (Hz)', 'Acumen Resampling Frequency (Hz)'}, ...
    'Analysis Settings', [1 50], ...
    {num2str(suggestedRef), num2str(suggestedEnd), '30', '25'});
if isempty(answer)
    error('Analysis settings dialog was cancelled.');
end
if exist("previewFig", "var") && isvalid(previewFig)
    close(previewFig);
end
cfg.refFrame         = round(str2double(answer{1}));
cfg.endFrame         = round(str2double(answer{2}));
cfg.dicFrameRateHz   = str2double(answer{3});
cfg.targetResampleHz = str2double(answer{4});

kin = computeCentroidKinematics(dic, cfg.insertionPointIDs, cfg.originPointIDs);

% Trim to the chosen analysis window; L0/displacement/strain are relative
% to the reference frame the user picked, not necessarily frame 1.
frameIdx     = cfg.refFrame:cfg.endFrame;
dicDistance     = kin.distance(frameIdx);
L0              = dicDistance(1);
dicDisplacement = dicDistance - L0;
dicStrain       = dicDisplacement / L0;
dicTime         = (0:numel(frameIdx)-1) / cfg.dicFrameRateHz;
fprintf("Initial DIC length L0 = %.4f\n", L0);

acumenStrain = acumen.AxialDisplacement / L0;   % same L0 normalization as the manual sheet

% Common time grid, clipped to where both signals actually overlap
tEnd  = min(dicTime(end), acumenTime(end));
tGrid = (0:1/cfg.targetResampleHz:tEnd)';

resampled.time         = tGrid;
resampled.dicDistance     = resampleSeries(dicTime, dicDistance, tGrid);
resampled.dicDisplacement = resampleSeries(dicTime, dicDisplacement, tGrid);
resampled.dicStrain       = resampleSeries(dicTime, dicStrain, tGrid);
resampled.acumenForce     = resampleSeries(acumenTime, acumen.AxialForce, tGrid);
resampled.acumenStrain    = resampleSeries(acumenTime, acumenStrain, tGrid);

plotResults(resampled);

resultsTable = table(resampled.time, resampled.dicDistance, resampled.dicDisplacement, ...
    resampled.dicStrain, resampled.acumenStrain, resampled.acumenForce, ...
    'VariableNames', {'Time_s', 'DIC_Distance', 'DIC_Displacement', 'DIC_Strain', 'Acumen_Strain', 'Acumen_Force_N'});
writetable(resultsTable, cfg.outputXlsx);

fprintf("results saved to dic_acumen_output folder\n");


%% ===================== LOCAL FUNCTIONS =====================

function dic = importDICRaw(filepath, sheet)
% Parses the "name | type | attribute | id | <frame1> | <frame2> | ..." DIC
% export layout. Each tracked point contributes three rows: pointx, pointy,
% pointz, with its id repeated in column 4. Accepts the raw .csv export as
% received from the DIC machine, or an .xlsx/.xls copy of the same data.

    [~, ~, ext] = fileparts(filepath);
    if ismember(lower(ext), [".xlsx", ".xls", ".xlsm"])
        raw = readcell(filepath, "Sheet", sheet);
    else
        raw = readcell(filepath, "FileType", "text", "Delimiter", ",");
    end
    attributeCol = raw(:, 3);
    idCol        = raw(:, 4);

    frameNumbers = cell2mat(raw(1, 5:end));
    numFrames    = numel(frameNumbers);

    isPointX = strcmp(attributeCol, "pointx");
    ids = unique(cell2mat(idCol(isPointX)));
    numPoints = numel(ids);

    X = nan(numPoints, numFrames);
    Y = nan(numPoints, numFrames);
    Z = nan(numPoints, numFrames);

    for k = 1:numel(ids)
        pid = ids(k);
        rowX = find(strcmp(attributeCol, "pointx") & cellfun(@(v) isequal(v, pid), idCol));
        rowY = find(strcmp(attributeCol, "pointy") & cellfun(@(v) isequal(v, pid), idCol));
        rowZ = find(strcmp(attributeCol, "pointz") & cellfun(@(v) isequal(v, pid), idCol));
        X(k, :) = cell2mat(raw(rowX, 5:end));
        Y(k, :) = cell2mat(raw(rowY, 5:end));
        Z(k, :) = cell2mat(raw(rowZ, 5:end));
    end

    dic.pointIDs     = ids;      % numPoints x 1, 0-indexed ids in row order of X/Y/Z
    dic.frameNumbers = frameNumbers;
    dic.X = X; dic.Y = Y; dic.Z = Z;   % each numPoints x numFrames
end

function kin = computeCentroidKinematics(dic, insertionIDs, originIDs)
% Averages the requested point IDs into two centroids per frame, then
% computes centroid-to-centroid distance for every frame. Reference frame,
% L0, displacement, strain, and the time vector all depend on the analysis
% window the user picks, so they're computed by the caller afterward.

    insertionRows = ismember(dic.pointIDs, insertionIDs);
    originRows    = ismember(dic.pointIDs, originIDs);

    cIns = [mean(dic.X(insertionRows, :), 1); mean(dic.Y(insertionRows, :), 1); mean(dic.Z(insertionRows, :), 1)];
    cOrg = [mean(dic.X(originRows, :), 1);    mean(dic.Y(originRows, :), 1);    mean(dic.Z(originRows, :), 1)];

    kin.centroidInsertion = cIns;  % 3 x numFrames
    kin.centroidOrigin    = cOrg;
    kin.distance = vecnorm(cIns - cOrg, 2, 1);          % 1 x numFrames
end

function daq = importAcumenDAQ(filepath)
% Reads an MTS TestSuite DAQ export. Handles two cases seen in practice:
%   - a genuinely delimited .txt file with a header block before the data
%   - a file saved with a .txt extension that is actually a zipped .xlsx
% Returns a struct with Time, AxialForce, AxialDisplacement (and CycleCount
% when present), all as column vectors in the file's original units.

    fid = fopen(filepath, "rb");
    firstBytes = fread(fid, 2, "uint8=>char")';
    fclose(fid);

    if isequal(firstBytes, "PK")
        tmpFile = [tempname, ".xlsx"];
        copyfile(filepath, tmpFile);
        raw = readcell(tmpFile);
        delete(tmpFile);
    else
        fid = fopen(filepath, "rt");
        lines = {};
        while true
            l = fgetl(fid);
            if ~ischar(l), break; end
            lines{end+1} = l; %#ok<AGROW>
        end
        fclose(fid);
        raw = cellfun(@(l) strsplit(l, "\t"), lines, "UniformOutput", false);
        maxCols = max(cellfun(@numel, raw));
        raw = cellfun(@(r) [r, repmat({''}, 1, maxCols - numel(r))], raw, "UniformOutput", false);
        raw = vertcat(raw{:});
    end

    % Find the header row: the first row containing a cell that starts with "Time"
    headerRowIdx = find(any(cellfun(@(c) ischar(c) && startsWith(strtrim(c), "Time"), raw), 2), 1);

    colNames = strtrim(string(raw(headerRowIdx, :)));
    dataStart = headerRowIdx + 2;   % skip the header row and the units row beneath it

    dataBlock = raw(dataStart:end, :);
    numericBlock = nan(size(dataBlock));
    for c = 1:size(dataBlock, 2)
        for r = 1:size(dataBlock, 1)
            v = dataBlock{r, c};
            if isnumeric(v)
                numericBlock(r, c) = v;
            elseif ischar(v) || isstring(v)
                numericBlock(r, c) = str2double(v);
            end
        end
    end
    validRows = ~all(isnan(numericBlock), 2);
    numericBlock = numericBlock(validRows, :);

    daq = struct();
    for c = 1:numel(colNames)
        name = colNames(c);
        if name == "" || ismissing(name), continue; end
        fieldName = regexprep(strtrim(name), "[^A-Za-z0-9]", "");
        daq.(fieldName) = numericBlock(:, c);
    end
end

function window = autoDetectRelaxationWindow(t, d, minHoldS, minDispFraction, settleTrimS, stillWindowS, stillFrac)
% Finds the longest constant-displacement hold at a non-trivial displacement
% level - i.e. the stress-relaxation segment, not the ramp/ cyclic-loading
% portions of the test. Same approach as FitCurve.m's autoDetectRelaxation.
% Returns [tStart tEnd] in the same time base as t, or [] if nothing qualifies.
    if nargin < 3, minHoldS = 15.0; end
    if nargin < 4, minDispFraction = 0.3; end
    if nargin < 5, settleTrimS = 0.3; end
    if nargin < 6, stillWindowS = 2.0; end
    if nargin < 7, stillFrac = 0.01; end

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

function yNew = resampleSeries(t, y, tNew)
% Linear interpolation of an irregularly/differently sampled series (t, y)
% onto a new common time grid tNew. Points outside the original range are NaN.
    [tUnique, iUnique] = unique(t);
    yNew = interp1(tUnique, y(iUnique), tNew, "linear");
end

function fig = plotResults(r)
    fig = figure("Name", "DIC + Acumen Strain/Force vs Time", "Color", "white");
    ax = axes(fig);
    hold(ax, "on");

    yyaxis(ax, "left");
    plot(ax, r.time, r.dicStrain, "-", "LineWidth", 1.2, "Color", [0.20 0.40 0.65], "DisplayName", "DIC Strain");
    plot(ax, r.time, r.acumenStrain, "-", "LineWidth", 1.2, "Color", [0.30 0.60 0.35], "DisplayName", "Acumen Strain");
    ylabel(ax, "Strain (mm/mm)");
    ax.YAxis(1).Color = [0.15 0.15 0.15];

    yyaxis(ax, "right");
    plot(ax, r.time, r.acumenForce, "-", "LineWidth", 1.2, "Color", [0.85 0.45 0.10], "DisplayName", "Axial Force");
    ylabel(ax, "Axial Force (N)");
    ax.YAxis(2).Color = [0.15 0.15 0.15];

    xlabel(ax, "Time (s)");
    title(ax, "Strain and Force vs Time");
    legend(ax, "show", "Location", "best", "Box", "off");

    ax.Box = "on";
    ax.GridColor = [0.85 0.85 0.85];
    ax.GridAlpha = 1;
    ax.LineWidth = 0.75;
    grid(ax, "on");
end
