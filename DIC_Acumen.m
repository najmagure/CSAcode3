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
cfg.dicFrameRateHz = 30;                    % DIC camera capture rate

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

cfg.targetResampleHz = 25;    % common frequency for the resampled/plotted data
cfg.outputXlsx       = fullfile(OUTPUT_DIR, "Resampled_Results.xlsx");

%% ===================== PIPELINE =====================
dic = importDICRaw(cfg.dicFile, cfg.dicSheet);

kin = computeCentroidKinematics(dic, cfg.insertionPointIDs, cfg.originPointIDs, cfg.dicFrameRateHz);
fprintf("Initial DIC length L0 = %.4f\n", kin.L0);

acumen = importAcumenDAQ(cfg.acumenFile);
if cfg.acumenTimeOffset == "auto"
    acumenTimeOffset = acumen.Time(1);
else
    acumenTimeOffset = cfg.acumenTimeOffset;
end
acumenTime   = acumen.Time - acumenTimeOffset;
acumenStrain = acumen.AxialDisplacement / kin.L0;   % same L0 normalization as the manual sheet

% Common time grid, clipped to where both signals actually overlap
tEnd  = min(kin.time(end), acumenTime(end));
tGrid = (0:1/cfg.targetResampleHz:tEnd)';

resampled.time         = tGrid;
resampled.dicDistance     = resampleSeries(kin.time, kin.distance, tGrid);
resampled.dicDisplacement = resampleSeries(kin.time, kin.displacement, tGrid);
resampled.dicStrain       = resampleSeries(kin.time, kin.strain, tGrid);
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

function kin = computeCentroidKinematics(dic, insertionIDs, originIDs, frameRateHz)
% Averages the requested point IDs into two centroids per frame, then
% computes centroid-to-centroid distance, displacement (relative to frame 1),
% and uniaxial engineering strain (displacement / initial distance).

    insertionRows = ismember(dic.pointIDs, insertionIDs);
    originRows    = ismember(dic.pointIDs, originIDs);

    cIns = [mean(dic.X(insertionRows, :), 1); mean(dic.Y(insertionRows, :), 1); mean(dic.Z(insertionRows, :), 1)];
    cOrg = [mean(dic.X(originRows, :), 1);    mean(dic.Y(originRows, :), 1);    mean(dic.Z(originRows, :), 1)];

    kin.centroidInsertion = cIns;  % 3 x numFrames
    kin.centroidOrigin    = cOrg;

    kin.distance = vecnorm(cIns - cOrg, 2, 1);          % 1 x numFrames
    kin.L0       = kin.distance(1);
    kin.displacement = kin.distance - kin.L0;
    kin.strain       = kin.displacement / kin.L0;

    numFrames = size(cIns, 2);
    kin.time = (0:numFrames-1) / frameRateHz;
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
