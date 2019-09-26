clear all
clc
clf
% hold on

%% Parameters:    
Period = 300 * 10^-9; % m
DPoints = 30; % k points along diagonal
sqPix = 11; %see below
Pix = sqPix; % # of mesh x
wmax = 1; % maximum plot value y-axis
eps = 8.85418782 * 10^-12; %F/m
DielectConst = 13; % Dielectirc constant of cylinder
rRatio = 0.39; % Ratio of radius of cylinder to period   
if rRatio > 0.4
    error('Settle down there cowboy.');
end
c =  2.99792458 * 10^8; %m/s
% eps = 8.85418782 * 10^-12; %F/m
Nmesh = Pix^2; % number of mesh elements in real space.

%% Create unit cell.  Use either CreatePost or CreateHole depending on the case
% CirclePixels = ones(sqPix);
%  CirclePixels = [1 1 1 1; 1 D_Const D_Const 1; 1 D_Const D_Const 1; 1 1 1 1];
[CirclePixels, radius] = CreateRegHexHole(Pix, DielectConst, rRatio);
%[CirclePixels, radius] = CreateRegHexPost(Pix, RefractiveIndex, rRatio);

%% Create epsilon matrix
epsA = EpsMatrix(CirclePixels); 

%% Create k points
% k = kPoints(Period, Period, DPoints); %(Rectangular path)
%k = kPointsReducedRegHexLattice(Period, DPoints); %(Hex Path)
k = kPointsReducedRegHexLatticeJBH(Period, DPoints); %(Hex Path)
k = cat(1, k(size(k, 1), :), k); % add gamma point to begining of k-vector
Nk = size(k, 1); % number of k vectors

% % plots k - path
% hold off;
% for i = 1:Nk
%     scatter(k(i, 1), k(i, 2));
%     pause(0.001);
%     hold on;
%     axis([min(k(:, 1)) max(k(:, 1)) min(k(:, 1)) max(k(:, 1))])
%     axis square
% end

%% Create k matrix
% FDFDmatrix = kMatrix(Period, Period, Pix, Pix, k);
FDFDmatrix = kMatrixRegHex(Period, Pix, k);

%% Solve for eigenvalues
eValues = cell(Nk, 1); % initializes matrix
eModes = cell(Nk, 1);
eModesNorm = cell(Nk, 1);
for m = 1:Nk
    eValues{m} = eig(FDFDmatrix{m}, epsA); % calculates Nmesh eigenvalues for each k point
    eModes{m} = sqrt(eValues{m}) * c;  % calculates the eModes = /omega (frequency)
    eModesNorm{m} = (eModes{m} * Period) / (2 * pi * c); % normalizes the frequency
end   
eModesNormMatrix = zeros(Nk, Nmesh); %I nstead of a Nk * 1 * Nmesh (3D) matrix, simplifing to a Nk * Nmesh (2D) matrix 
for j = 1:Nk
    eModesNormMatrix(j, :) = eModesNorm{j}'; % makes the 2D matrix
end

%% plot band structure
kPlotVector = zeros(Nk, 1); % initializes x-axis (k)
for i = 1:Nk
    kPlotVector(i) = i;
end
%   eModesNormMatrix = cat(2, kPlotVector, eModesNormMatrix); 
%   PlotModes = plot(eModesNormMatrix(:, 1), eModesNormMatrix(:, 2:(Pix * Pix + 1)));
%   PlotModes = plot(eModesNormMatrix(:, 1), eModesNormMatrix(:, 2:Nmesh));
PlotModes = zeros(Nk, 1);
for i = 1:(Nmesh) 
    PlotModes(i) = plot(kPlotVector, eModesNormMatrix(:, i),'.');
    hold on;
end

% create xaxis labels
Mpoint = 1 + round( sqrt(3) * DPoints / 2 );
Kpoint = round ( Mpoint + DPoints / 2 );
set(gca, 'Xtick', [1, Mpoint, Kpoint, Nk]); 
set(gca, 'XTickLabel', ['G'; 'M'; 'K'; 'G']);
% %     text(1, -0.075 * 10^16, '\Gamma');
% %     text(Nk, -0.075 * 10^16, '\Gamma');
axis([1 Nk 0 1])
hold on 

% draw vertical lines
plot([Mpoint Mpoint], [0 wmax], 'k');
plot([Kpoint Kpoint], [0 wmax], 'k');
disp('done');