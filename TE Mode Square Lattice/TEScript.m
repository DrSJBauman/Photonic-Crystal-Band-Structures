%The output of this function will be the A matrix for the eigenvalue problem (FDFDmatrix)
%It is currently coded as a script for ease of testing.
clear
clc
clf

%% Inputs
xPeriod = 300 * 10^-9; % Period between atoms along basis vector a1 [m]
yPeriod = 300 * 10^-9; % Period between atoms along basis vector a2 [m]
DPoints = 50; % Number of k points along the diagonal in k-space
sqPix = 5; % # of mesh points in one unit cell of the PC in one direction [elements]
xPix = sqPix; % # of mesh points in x [elements]
yPix = sqPix; % Sloppy coding [elements]
Pix = xPix; %More sloppy/goofy coding [elements]
wmax = 1; % Maximum plot value for the y-axis (Normalized frequency)
EpsC = 8.9; % Relative dielectric constant of post/hole [unitless]
EpsB = 1; % Relative dielectric constant of bulk [unitless]
theta = 90; % Angle between basis vectors [degrees]
rRatio = 0.4; % Ratio of radius of cylinder to period
ax = xPeriod / xPix; % Size of the x mesh elements [m]
ay = yPeriod / yPix; % Size of the y mesh elements [m]

%% Setting up the dielectric matrix so we can access each mesh element's dielectric value
%Creating a matrix to contain the dielectric constant values for each mesh point in the unit cell 
%instead of storing the coordinates and dielectric values in a (n * 3) matrix.
UnitCell = GenGeom(xPeriod, yPeriod, theta, xPix, yPix, EpsC, EpsB, rRatio); %Calling the geometry function
k = kPointsReduced(xPeriod, yPeriod, DPoints); %Calling the k points function

% make UnitCell into Matrix
UnitRow = 1 ./ UnitCell(:, 3);
p = zeros(Pix); %Setting up a matrix to contain the dielectric values at positions corresponding to their mesh positions in the unit cell
for i = 1:Pix
    p(i, :) = UnitRow(((i - 1) * Pix + 1):(i * Pix));
end

% for i = 1:Pix
%     for j = 1:Pix
%         p(i,j)=10*rand();
%     end
% end

H_Top = cell(numel(k(:, 1)) , 1);
H_Bot = cell(numel(k(:, 1)) , 1);
H = cell(numel(k, 1), 1);  

%% Looping through every k point
for m = 1:numel(k(:, 1)) 
    %% Hxi (Top left)
    % Initializes the (Pix^2 * Pix^2) matrix that contains Pix number of (Pix * Pix) matrices
    Hxi = cell(Pix);
    for i = 1:Pix
        for j = 1:Pix
            Hxi{i, j} = zeros(Pix);
        end
    end
    
     for i = 1:Pix % Loop for diagonal of Pix * (Pix * Pix) Matrices
        % Main diagonal
        for id = 1:Pix
           Hxi{i, i}(id, id) = p(i, id) * (2 / ax^2 + 2 / ay^2);
        end

        % Upper diagonal
        Hxi{i, i}(1, 2) = (-p(i, 2) + p(i, Pix) - 4 * p(i, 1)) / (4 * ay^2); 
        for id = 2:(Pix - 1)
           Hxi{i, i}(id, id + 1) = (-p(i, id + 1) + p(i, id - 1) - 4 * p(i, id)) / (4 * ay^2); 
        end

        % Lower diagonal
        Hxi{i, i}(2, 1) = (p(i, Pix) - p(i, 1) - 4 * p(i, 2)) / (4 * ay^2);
        for id = 2:(Pix - 1)
           Hxi{i, i}(id + 1, id) = (p(i, id - 1) - p(i, id) - 4 * p(i, id + 1)) / (4 * ay^2);
        end
        
        % Top corner
        Hxi{i, i}(1, Pix) = exp(1i * k(m, 2) * yPeriod) * (p(i, 2) - p(i, Pix) - 4 * p(i, 1)) / (4 * ay^2);

        % Bot corner
        Hxi{i, i}(Pix, 1) = exp(-1i * k(m, 2) * yPeriod) * (-p(i, 1) + p(i, 2) - 4 * p(i, Pix)) / (4 * ay^2);
     end
     
     % Upper and lower diagonal of (Pix-1) by (Pix * Pix) Matrices
     for i = 1:(Pix - 1)
         for id = 1:Pix
            Hxi{i, i + 1}(id, id) = -p(i, id) / ax^2;
            Hxi{i + 1, i}(id, id) = -p(i + 1, id) / ax^2;
         end
     end
     
     % Corner Pix * Pix matrix
     for id = 1:Pix
         Hxi{1, Pix}(id, id) = exp(1i * k(m, 1) * xPeriod) * (-p(1, id) / ax^2);
         Hxi{Pix, 1}(id, id) = exp(-1i * k(m, 1) * xPeriod) * (-p(Pix, id) / ax^2);
     end
     
     
    %% Hyj (Bottom Right)
    % Initializes the (Pix^2 * Pix^2) matrix that contains Pix number of (Pix * Pix) matrices
    Hyj = cell(Pix);
    for i = 1:Pix
        for j = 1:Pix
            Hyj{i, j} = zeros(Pix);
        end
    end
     for i = 1:Pix % Loop for diagonal of Pix by (Pix * Pix) Matrices

        % Main diagonal
        for id = 1:Pix
           Hyj{i, i}(id, id) = p(i, id) * (2 / ax^2 + 2 / ay^2);
        end

        % Upper diagonal
        for id = 1:(Pix - 1)
           Hyj{i, i}(id, id + 1) = -p(i, id) / ay^2;
        end

        % Lower diagonal
        for id = 2:Pix
           Hyj{i, i}(id, id - 1) = -p(i, id) / ay^2;
        end
        
        % Top corner
         Hyj{i, i}(1, Pix) = exp(1i * k(m, 2) * yPeriod) * (-p(i, 1) / ay^2);

        % Bot corner
         Hyj{i, i}(Pix, 1) = exp(-1i * k(m, 2) * yPeriod) * (-p(i, Pix) / ay^2);
     end
     
     % Upper and lower diagonal of (Pix-1) by (Pix * Pix) Matrices
     for id = 1:Pix
            Hyj{1, 2}(id, id) = (-p(2, id) + p(Pix, id) - 4 * p(1, id)) / (4 * ax^2);
            Hyj{2, 1}(id, id) = (p(Pix, id) - p(1, id) - 4 * p(2, id)) / (4 * ax^2);
     end
         
     for i = 2:(Pix - 1)
         for id = 1:Pix
            Hyj{i, i + 1}(id, id) = (-p(i + 1, id) + p(i - 1, id) - 4 * p(i, id)) / (4 * ax^2); 
            Hyj{i + 1, i}(id, id) = (p(i - 1, id) - p(i, id) - 4 * p(i + 1, id)) / (4 * ax^2);
         end
     end
     
     % Corner Pix * Pix matrix
     for id = 1:Pix
         Hyj{1, Pix}(id, id) = exp(1i * k(m, 1) * xPeriod) * (p(2, id) - p(Pix, id) - 4 * p(1, id)) / (4 * ax^2);
         Hyj{Pix, 1}(id, id) = exp(-1i * k(m, 1) * xPeriod) * (p(1, id) + p(2, id) - 4 * p(Pix, id)) / (4 * ax^2);
     end
     
    %% Hyi (Top Right)
    % Initializes the (Pix^2 * Pix^2) matrix that contains Pix number of (Pix * Pix) matrices
    Hyi = cell(Pix);
    for i = 1:Pix
        for j = 1:Pix
            Hyi{i, j} = zeros(Pix);
        end
    end
    
    % Upper and lower diagonal of (Pix-1) by (Pix * Pix) matrices in size
    for i = 1:(Pix - 1)
        HyiL{i, i + 1}(Pix, Pix) = p(i, 1);
        HyiR{i, i + 1}(1, 1) = p(i, Pix);
        for id = 1:(Pix - 1)
            HyiL{i, i + 1}(id, id) = p(i, id + 1);
            HyiR{i, i + 1}(id + 1, id + 1) = p(i, id);
        end
        Hyi{i, i + 1} = (HyiL{i, i + 1} - HyiR{i, i + 1}) / (4 * ax * ay);
    end
    
    % Bottom diagonal
    for i = 1:(Pix - 1)
        HyiL{i + 1, i}(Pix, Pix) = p(i + 1, 1);
        HyiR{i + 1, i}(1, 1) = p(i + 1, Pix);
        for id = 1:(Pix - 1)
            HyiL{i + 1, i}(id, id) = p(i + 1, id + 1);
            HyiR{i + 1, i}(id + 1, id + 1) = p(i + 1, id);
        end
        Hyi{i + 1, i} = (-HyiL{i + 1, i} + HyiR{i + 1, i}) / (4 * ax * ay);
    end
    
    % Top right Pix * Pix matrix
    HyiL{1, Pix}(Pix, Pix) = p(1, 1);
    HyiR{1, Pix}(1, 1) = p(1, Pix);
    for id = 1:(Pix - 1)
        HyiL{1, Pix}(id, id) = p(1, id + 1);
        HyiR{1, Pix}(id + 1, id + 1) = p(1, id);
    end
    Hyi{1, Pix} = exp(1i * k(m, 1) * xPeriod) * (-HyiL{1, Pix} + HyiR{1, Pix}) / (4 * ax * ay);
    
    % Bottom Left Pix * Pix matrix
    HyiL{Pix, 1}(Pix, Pix) = p(Pix, 1);
    HyiR{Pix, 1}(1, 1) = p(Pix, Pix);
    for id = 1:(Pix - 1)
        HyiL{Pix, 1}(id, id) = p(Pix, id + 1);
        HyiR{Pix, 1}(id + 1, id + 1) = p(Pix, id);
    end
    Hyi{Pix, 1} = exp(-1i * k(m, 1) * xPeriod) * (HyiL{Pix, 1} - HyiR{Pix, 1}) / (4 * ax * ay);
   
    
    %% Hxj (Bottom Left)
    % Initializes the (Pix^2 * Pix^2) matrix that contains Pix number of (Pix * Pix) matrices
    Hxj = cell(Pix);
    for i = 1:Pix
        for j = 1:Pix
            Hxj{i, j} = zeros(Pix); % Initializes the (Pix * Pix) matrix
        end
    end
    
    % Top Right and Bottom Left corners of Pix by (Pix * Pix) matrices
    HxjL{Pix, Pix}(1, Pix) = p(1, 1);  
    HxjR{1, 1}(1, Pix) = p(Pix, 1);
    HxjL{Pix, Pix}(Pix, 1) = p(1, Pix);  
    HxjR{1, 1}(Pix, 1) = p(Pix, Pix);
    for i = 1:(Pix - 1)
        HxjL{i, i}(1, Pix) = p(i + 1, 1);
        HxjR{i + 1, i + 1}(1, Pix) = p(i, 1);    
        HxjL{i, i}(Pix, 1) = p(i + 1, Pix);
        HxjR{i + 1, i + 1}(Pix, 1) = p(i, Pix);    
    end
    for i = 1:Pix
        Hxj{i, i}(1, Pix) = exp(1i * k(m, 2) * yPeriod) * (-HxjL{i, i}(1, Pix) + HxjR{i, i}(1, Pix)) / (4 * ax * ay);
        Hxj{i, i}(Pix, 1) = exp(-1i * k(m, 2) * yPeriod) * (HxjL{i, i}(Pix, 1) - HxjR{i, i}(Pix, 1)) / (4 * ax * ay);
    end
    
    % Upper and bottm diagonals of Pix by (Pix * Pix) matrices
    % THIS IS WHERE WE LEFT OFF
    %Bauman - None of these can be left outside of the loops because either id or i is always looping in all cases
    %Check above, where we made p a random integer, because I turned it off while cleaning up the code today.
    for i = 1:(Pix - 1)
        for id = 1:(Pix - 1)
            HxjL{i, i}(id, id + 1) = p(i + 1, id);
            HxjR{i + 1, i + 1}(id, id + 1) = p(i, id);
            HxjL{i, i}(id + 1, id) = p(i + 1, id + 1);
            HxjR{i + 1, i + 1}(id + 1, id) = p(i + 1, id + 1);
            HxjL{Pix, Pix}(id, id + 1) = p(1, id);
            HxjR{1, 1}(id, id + 1) = p(Pix, id);
            HxjL{Pix, Pix}(id + 1, id) = p(1, id+1);
            HxjR{1, 1}(id + 1, id) = p(Pix, id+1);
        end
    end
    for i = 1:Pix
        for id = 1:(Pix - 1)
            Hxj{i, i}(id, id + 1) = (HxjL{i, i}(id, id + 1) - HxjR{i, i}(id, id + 1)) / (4 * ax * ay);
            Hxj{i, i}(id + 1, id) = (-HxjL{i, i}(id + 1, id) + HxjR{i, i}(id + 1, id)) / (4 * ax * ay);
        end
    end
    %STOP HERE
    
    %% Concatenating the 4 (Pix^2 * Pix^2) matrices into a 2 * (Pix^2 by Pix^2) matrix
    Hxi = cell2mat(Hxi);
    Hyi = cell2mat(Hyi);
    Hxj = cell2mat(Hxj);
    Hyj = cell2mat(Hyj);
    H_Top{m} = cat(2, Hxi, Hyi);
    H_Bot{m} = cat(2, Hxj, Hyj);
    H{m} = cat(1, H_Top{m}, H_Bot{m});
end

FDFDmatrix = H'; % Transpose

spy(cell2mat(Hxj))