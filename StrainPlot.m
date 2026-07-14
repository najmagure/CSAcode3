clear
clc
close all

%Select DIC File to upload
[fileDIC,pathDIC] = uigetfile({'*.csv;*.xlsx'},...
    'Select DIC File');

if isequal(fileDIC,0)
    error('No DIC file selected');
end

dicFile = fullfile(pathDIC,fileDIC);

[~,~,ext] = fileparts(dicFile);

if strcmpi(ext,'.csv')
    raw = readcell(dicFile);
else
    raw = readcell(dicFile);
end

DIC = cell2mat(raw(3:end,5:end));

%Insertion Centroid Calculation

xIns = mean(DIC([1 4 7],:),1);
yIns = mean(DIC([2 5 8],:),1);
zIns = mean(DIC([3 6 9],:),1);

%Origin Centroid Calculation

xOrg = mean(DIC([19 22 25],:),1);
yOrg = mean(DIC([20 23 26],:),1);
zOrg = mean(DIC([21 24 27],:),1);

%Distance Calculation

distance = sqrt( ...
    (xIns - xOrg).^2 + ...
    (yIns - yOrg).^2 + ...
    (zIns - zOrg).^2 );
%User Input

answer = inputdlg( ...
    {'Reference Frame (T = 0)', ...
    'Last Frame to Analyze', ...
    'DIC Frequency (Hz)', ...
    'Acumen Resampling Frequency (Hz)'}, ...
    'Analysis Settings', ...
    [1 50]);

refFrame = str2double(answer{1});
endFrame = str2double(answer{2});
dicFs = str2double(answer{3});
targetFs = str2double(answer{4});

%Setting range for DIC
distance = distance(refFrame:endFrame);

%Initial Distance L0

L0 = distance(1);

fprintf('\nReference Length (L0) = %.4f mm\n',L0);
%DIC Displacement & Strain

dicDisplacement = distance - L0;

dicStrain = dicDisplacement ./ L0;

%DIC Time Vector


nFrames = length(distance);

dicTime = (0:nFrames-1)'/dicFs;


%Select Acumen File to upload
[fileAC,pathAC] = uigetfile({'*.csv;*.xlsx'},...
    'Select Acumen File');

if isequal(fileAC,0)
    error('No Acumen file selected');
end

ACFile = fullfile(pathAC,fileAC);

%Read Acumen Data


ACTable = readtable(ACFile);

disp('Available Acumen Columns:')
disp(ACTable.Properties.VariableNames)

% Modify if needed

acDisp = ACTable.AxialDisplacement;
acForce = ACTable.AxialForce;


%Resampling

originalFs = 100;

[p,q] = rat(targetFs/originalFs);

dispResampled  = resample(acDisp,p,q);
forceResampled = resample(acForce,p,q);

%Time Vector

timeResampled = (0:length(dispResampled)-1)'/targetFs;

%Calculate Acumen strain
acumenStrain = dispResampled ./ L0;


figure

plot(dicTime,dicStrain,...
    'LineWidth',2)

hold on

plot(timeResampled,acumenStrain,...
    'LineWidth',2)

xlabel('Time (s)')
ylabel('Strain')
title('DIC Strain vs Acumen Strain')
legend('DIC','Acumen')
grid on

%%Currently stretches/compresses DIC Data to fit the number of frames
%%into set test duration, need to have it either assign time based on set 0
%%and then align frames to frequency(prob best) or just change it in the
%%file.
%Need to do the same to Acumen it currently works using time data from
%acumen file
%Need to also be able to set intial distance based on the frame/time used
%to differentiate each test but might not have to code for it if I choose
%the starting frame from the beginning