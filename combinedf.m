%% 1. Load Specimen
[file, path] = uigetfile({'*.stl;*.STL','STL files (*.stl)'; '*.*','All files'}, 'Select STL file');

fname = fullfile(path,file);

TR = stlread(fname);
V = TR.Points;
F = TR.ConnectivityList;

%% 2. PCA Main Axis 
[coeff, ~, ~, ~, ~, mu] = pca(V);

n = coeff(:,1)' / norm(coeff(:,1));

%% 3. Cut Positions
perc = [20, 50, 80];

proj = V * n'; 

lengthAlongN = max(proj) - min(proj);
minp = min(proj);
p_coords = minp + (perc(:)/100) * lengthAlongN;

P0_proj = mu * n';

%% 4. Cross-Section Calculation
num = numel(p_coords);
areas = zeros(num,1);
cutPoints = zeros(num,3);
loopsByCut = cell(num,1);

for k = 1:num
    Ppos = mu + (p_coords(k) - mu*n') * n; 
    cutPoints(k,:) = Ppos;

    [areas(k), loopsByCut{k}, ~] = crossSectionAreaFromMesh(V, F, Ppos, n);
end

%% 5. Display Results
T = table(perc(:), p_coords, areas, ...
    'VariableNames', {'Percent', 'ProjCoord', 'Area'});
disp(T);

%% 6. Plot Setup
figure;
trisurf(F, V(:,1), V(:,2), V(:,3), ...
    'FaceAlpha', 0.6, ...
    'EdgeColor', 'none', ...
    'DisplayName', 'Specimen');
hold on;
axis equal;
colors = lines(num);

%% 7. In-Plane Directions
if abs(n(3)) < 0.9
    tmp = [0 0 1];
else
    tmp = [0 1 0];
end

u = cross(n,tmp); u = u / norm(u);
v = cross(n,u); v = v / norm(v);

scale = 1.2 * max(range(V));

%% 8. PCA Axis Visualization
nLine = [mu - 0.55*lengthAlongN*n; mu + 0.55*lengthAlongN*n];

plot3(nLine(:,1), nLine(:,2), nLine(:,3), 'k--', ...
    'LineWidth', 2, ...
    'DisplayName', 'n axis');

%% 9. Draw Cuts
for k = 1:num
    Ppos = cutPoints(k,:);
    loops3D = loopsByCut{k};

    for L = 1:numel(loops3D)
        pts = loops3D{L};

        if L == 1
            plot3(pts(:,1), pts(:,2), pts(:,3), '-', ...
                'Color', colors(k,:), ...
                'LineWidth', 1.5, ...
                'DisplayName', sprintf('%d%% cut', perc(k)));
        else
            plot3(pts(:,1), pts(:,2), pts(:,3), '-', ...
                'Color', colors(k,:), ...
                'LineWidth', 1.5, ...
                'HandleVisibility', 'off');
        end
    end

    [su,sv] = meshgrid([-1 1], [-1 1]);
    planeCorners = Ppos + (su(:)*scale).*u + (sv(:)*scale).*v;

    patch('Vertices', planeCorners, ...
        'Faces', [1 2 4 3], ...
        'FaceColor', colors(k,:), ...
        'FaceAlpha', 0.15, ...
        'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
end

%% 10. Final Format
title('Specimen Cross Sections from PCA Axis');
xlabel('X'); ylabel('Y'); zlabel('Z');
legend('show', 'Location', 'best');
grid on;
view(3);
camlight;
lighting gouraud;

function [areaTotal, loops3D, areas] = crossSectionAreaFromMesh(V, F, P0, n)

n = n(:)' / norm(n);
tol = 1e-12;
edges = [1 2; 2 3; 3 1];

% collect intersection segments
segments = zeros(size(F,1)*2,6);
segcount = 0;
for i = 1:size(F,1)
    tri = V(F(i,:),:);                 % 3x3
    d = (tri - P0) * n';               % signed distances
    pts = zeros(0,3);
    for e = 1:3
        a = edges(e,1); b = edges(e,2);
        da = d(a); db = d(b);
        if abs(da) < tol && abs(db) < tol
            pts = [pts; tri(a,:); tri(b,:)]; 
        elseif abs(da) < tol
            pts = [pts; tri(a,:)]; 
        elseif abs(db) < tol
            pts = [pts; tri(b,:)]; 
        elseif da * db < 0
            t = da / (da - db);
            pnt = tri(a,:) + t*(tri(b,:) - tri(a,:));
            pts = [pts; pnt]; %#ok<AGROW>
        end
    end
    if size(pts,1) >= 2
        pts = unique(round(pts,12),'rows');
        if size(pts,1) == 2
            segcount = segcount + 1;
            segments(segcount,:) = [pts(1,:) pts(2,:)];
        elseif size(pts,1) > 2
            for k = 1:size(pts,1)-1
                segcount = segcount + 1;
                segments(segcount,:) = [pts(k,:) pts(k+1,:)];
            end
        end
    end
end
segments = segments(1:segcount,:);

% no intersection
if isempty(segments)
    areaTotal = 0;
    loops3D = {};
    areas = [];
    return;
end

% unique points and adjacency
P_list = [segments(:,1:3); segments(:,4:6)];
P_list = unique(round(P_list,12),'rows','stable');

toIndex = @(pnt) find(all(abs(P_list - pnt) < 1e-8, 2), 1);
E = zeros(size(segments,1),2);
for k = 1:size(segments,1)
    p1 = segments(k,1:3); p2 = segments(k,4:6);
    E(k,1) = toIndex(p1); E(k,2) = toIndex(p2);
end

% stitch segments into loops
used = false(size(E,1),1);
loopsIdx = {};
while any(~used)
    e = find(~used,1); used(e) = true;
    chain = E(e,:);
    startIdx = chain(1); endIdx = chain(2);
    extended = true;
    while extended
        extended = false;
        for ee = find(~used)'
            if E(ee,1) == endIdx
                used(ee) = true; endIdx = E(ee,2); chain = [chain, endIdx]; extended = true; break;
            elseif E(ee,2) == endIdx
                used(ee) = true; endIdx = E(ee,1); chain = [chain, endIdx]; extended = true; break;
            end
        end
        if ~extended
            for ee = find(~used)'
                if E(ee,2) == startIdx
                    used(ee) = true; startIdx = E(ee,1); chain = [startIdx, chain]; extended = true; break;
                elseif E(ee,1) == startIdx
                    used(ee) = true; startIdx = E(ee,2); chain = [startIdx, chain]; extended = true; break;
                end
            end
        end
    end
    if chain(1) ~= chain(end)
        chain = [chain, chain(1)];
    end
    loopsIdx{end+1} = chain;
end

% compute areas by projecting into plane basis
% construct in-plane orthonormal basis u,v
u_guess = P_list(1,:) - P0;
if norm(u_guess) < 1e-8
    u_guess = P_list(min(2,size(P_list,1)),:) - P0;
end
u = u_guess - n*(n*u_guess');
u = u / norm(u);
v = cross(n,u); v = v / norm(v);

numLoops = numel(loopsIdx);
loops3D = cell(numLoops,1);
areas = zeros(numLoops,1);
for k = 1:numLoops
    idxs = loopsIdx{k};
    pts3 = P_list(idxs,:);
    if isequal(pts3(1,:), pts3(end,:)), pts3(end,:) = []; end
    loops3D{k} = pts3;
    XY = [(pts3 - P0) * u', (pts3 - P0) * v'];
    c = mean(XY,1);
    ang = atan2(XY(:,2)-c(2), XY(:,1)-c(1));
    [~, order] = sort(ang);
    XYs = XY(order,:);
    areas(k) = polyarea(XYs(:,1), XYs(:,2));
end
areaTotal = sum(areas);
end