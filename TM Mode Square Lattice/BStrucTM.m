clear all
clc
clf

%% Parameters:   
Period = 275 * 10^-9; % Distance between atoms [m]
P1 = Period; % Distance along basis vector 1 [m]
P2 = Period; % Distance along basis vector 2 [m]
DPoints = 50; % Number of k points along the diagonal of the calculated irreducible BZ
Pix = 50; % Number of mesh points in the photonic crystal unit cell in one direction
xPix = Pix; % Number of mesh points along x
yPix = Pix; % Number of mesh points along y
wmax = 0.8; % Maximum y-axis value on the final band structure plot
EpsC = 11.8; % Relative dielectric constant of the Cylinders (posts/holes) [unitless]
EpsB = 1; % Relative dielectric constant of the Bulk surrounding the posts/holes [unitless]
theta = 90; % Angle between the real space basis vectors [degrees]
rRatio = 0.25; % Ratio between the radius of the cylinders and the period   
if rRatio > 0.49
    error('Settle down there cowboy. Try a smaller rRatio.');
end
c =  2.99792458 * 10^8; % Speed of light [m/s]
Nmesh = xPix * yPix; % Number of total mesh elements in the unit cell in real space

%% Create unit cell by applying dielectric constant values to mesh positions.
UnitCell = GenGeom(P1, P2, theta, xPix, yPix, EpsC, EpsB, rRatio);
%Use either CreatePost or CreateHole depending on the case
% [UnitCell, radius] = CreateHole(xPix, yPix, EpsC, rRatio);
% [CirclePixels, radius] = CreatePost(xPix, yPix, EpsC, rRatio);

%% Create epsilon matrix
% epsA = EpsMatrix(UnitCell); 
epsA = zeros(Pix^2);
for i = 1:Pix^2
    epsA(i, i) = UnitCell(i, 3); %Creates a diagonal matrix of epsilon values.
end

%% Create k points
k = kPointsReduced(P1, P2, DPoints); % Square path in k-space
k = cat(1, k(size(k, 1), :), k); % Add gamma point to begining of k-vector
Nk = size(k, 1); % Number of k vectors

% plots k path
% hold off;
% for i = 1:Nk
%     scatter(k(i, 1), k(i, 2));
%     pause(0.001);
%     hold on;
%     axis([min(k(:, 1)) max(k(:, 1)) min(k(:, 2)) max(k(:, 2))])
% end

%% Create cell containing the A matrices
% FDFDmatrix is an array of Nk number of A matrices (in the generalized eigenvalue problem)
% There is a separate one for each k vector.
FDFDmatrix = kMatrix(P1, P2, xPix, yPix, k);
% See notes above and https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix#Generalized_eigenvalue_problem

%% Solve for eigenvalues
eValues = cell(Nk, 1); % Initializing matrices
eModes = cell(Nk, 1);
eModesNorm = cell(Nk, 1);
for m = 1:Nk
    eValues{m} = eig(FDFDmatrix{m}, epsA); % Calculates Nmesh eigenvalues for each k point
    eModes{m} = sqrt(eValues{m}) * c;  % Calculates the eModes = omega (angular frequency)
    eModesNorm{m} = (eModes{m} * P1) / (2 * pi * c); % Normalizes the frequency
    lbdaModes{m} = (Period * 1e9 ./ eModesNorm{m}); % Solve for the wavelengths to use in COMSOL comparisons
end   
eModesNormMatrix = zeros(Nk, Nmesh); %Instead of a Nk * 1 * Nmesh (3D) matrix, simplifying to a Nk * Nmesh (2D) matrix to contain normalized freq values
lbda = zeros(Nk, Nmesh); % Instead of a Nk * 1 * Nmesh (3D) matrix, making a Nk * Nmesh (2D) matrix to contain wavelength values
for j = 1:Nk % Cycling through all k points
    eModesNormMatrix(j, :) = eModesNorm{j}'; % Matrix with a row of normalized frequency values for each k point
    lbda(j, :) = lbdaModes{j}'; % Matrix with a row of wavelength values for each k point
end

%% Plot band structure
kPlotVector = zeros(Nk, 1); % Initializes x-axis (k)
for i = 1:Nk
    kPlotVector(i) = i;
end

PlotModes = zeros(Nk, 1);
for i = 1:Nmesh % Cycle through all of the solved eigenvalues for a given k
    PlotModes(i) = plot(kPlotVector, eModesNormMatrix(:, i), 'b'); %Normalized
    hold on;
end

% Create x-axis labels for high symmetry k-points
Xpoint = round(DPoints / sqrt(2));
Spoint = 2 .* Xpoint;
set(gca, 'Xtick', [1, Xpoint, Spoint, Nk]); 
set(gca, 'XTickLabel', ['G'; 'X'; 'S'; 'G']);
axis([1 Nk 0 wmax])
hold on 

% Make the vertical lines for high symmetry k-points
plot([Xpoint Xpoint], [0 wmax], 'k');
plot([Spoint Spoint], [0 wmax], 'k');
hold off

%% Plot band structure and COMSOL results together for comparison
figure
% lbda = real(2 * pi * (c ./ eModesNormMatrix) * Period); % Solve for the wavelengths to use in COMSOL comparisons
if EpsB == 1
    ComsolComparePosts(Nk, Nmesh, lbda, Xpoint);
elseif EpsC == 1
    ComsolCompareHoles(Nk, Nmesh, lbda, Xpoint);
end

%% Plot to test the shift in wavelength versus mesh size
% figure
% GapShiftTest(Nk, Pix, Nmesh, lbda, Xpoint, Spoint);

%% Use this to export data to FDFD Pretty Plots Excel sheet if Nk=121
ExcelData = zeros(size(lbda(:,1)),6);
for i = 1:6
    ExcelData(:,i) = real(lbda(:,i)); %After running, open the variable ExcelData and copy/paste values into Excel
end

%%
disp('Code has finished running.');