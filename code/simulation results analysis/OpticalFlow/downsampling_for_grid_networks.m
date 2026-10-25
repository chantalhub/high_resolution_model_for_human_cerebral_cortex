% https://github.com/BrainDynamicsUSYD/BrainVortexToolbox/blob/main/Sub-function/Preprocessing_main.m

% left hemisphere
% step 1
position = double(loc);
k = alphaShape(position(1:29556, 1), position(1:29556, 2));
[a, b] = k.boundaryFacets();
xCord = -250: 2: 250; 
yCord = -150: 2: 200;
bw = poly2mask(b(:, 1) - min(xCord) + 1, b(:, 2) - min(yCord) + 1, max(yCord) - min(yCord) + 1, max(xCord) - min(xCord) + 1);
area = bw(211:233, 54:69);

% step 2
position = double(loc);
k = alphaShape(position(1:29556, 1), position(1:29556, 2), 8.2);
[a, b] = k.boundaryFacets();
xCord = -250: 2: 250; 
yCord = -150: 2: 200;
bw = poly2mask(b(:, 1) - min(xCord) + 1, b(:, 2) - min(yCord) + 1, max(yCord) - min(yCord) + 1, max(xCord) - min(xCord) + 1);

bw_modified = bw;
bw_modified(211:233, 54:69) = area;

mask = double(bw_modified(1:2:end, 1:2:end));
mask(mask == 0) = nan;
[X, Y] = meshgrid(xCord, yCord);
sigSSpec_reshape = [];
tic
if sum(mask(:)) ~= 0
    flagMask = 1;
end
for iTime = 1:size(firing_rate, 2)
    sigSSpec_reshape(:, :, iTime) = griddata(double(loc(1:29556, 1)), double(loc(1:29556, 2)), firing_rate(1:29556, iTime), X, Y, 'cubic');
    warning off;
    if flagMask
        sigSSpec_reshape(:, :, iTime) = sigSSpec_reshape(:, :, iTime) .* mask;
        warning off;
    end
end
toc


%{
% right hemisphere
position = double(loc);
k = alphaShape(position(29557:end, 1), position(29557:end, 2), 4);
[a, b] = k.boundaryFacets();
xCord = -270: 2: 230; 
yCord = -180: 2: 170;
bw = poly2mask(b(:, 1) - min(xCord) + 1, b(:, 2) - min(yCord) + 1, max(yCord) - min(yCord) + 1, max(xCord) - min(xCord) + 1);
area = bw(251:274, 444:455);

% step 2
position = double(loc);
k = alphaShape(position(29557:end, 1), position(29557:end, 2), 8.2);
[a, b] = k.boundaryFacets();
xCord = -270: 2: 230; 
yCord = -180: 2: 170;
bw = poly2mask(b(:, 1) - min(xCord) + 1, b(:, 2) - min(yCord) + 1, max(yCord) - min(yCord) + 1, max(xCord) - min(xCord) + 1);

bw_modified = bw;
bw_modified(251:274, 444:455) = area;

mask = double(bw_modified(1:2:end, 1:2:end));
mask(mask == 0) = nan;
[X, Y] = meshgrid(xCord, yCord);
sigSSpec_reshape = [];
tic
if sum(mask(:)) ~= 0
    flagMask = 1;
end
for iTime = 1:size(firing_rate, 2)
    sigSSpec_reshape(:, :, iTime) = griddata(double(loc(29557:end, 1)), double(loc(29557:end, 2)), firing_rate(29557:end, iTime), X, Y, 'cubic');
    warning off;
    if flagMask
        sigSSpec_reshape(:, :, iTime) = sigSSpec_reshape(:, :, iTime) .* mask;
        warning off;
    end
end
toc
%}
