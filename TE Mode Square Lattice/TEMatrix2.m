function [FDFDmatrix, p] = TEMatrix2(xPeriod, yPeriod, xPix, yPix, k, UnitCell)
% Creates A (FDFDmatrix) = output
% xPeriod = unit cell period in x
% yPeriod = unit cell period in y
% xPix = number of mesh points in x direction
% yPix = number of mesh points in y direction
% UnitCell = Matrix of each mesh point of unit cell
Pix = xPix; %More sloppy/goofy coding [elements]
ax = xPeriod / xPix; % Size of the x mesh elements [m]
ay = yPeriod / yPix; % Size of the y mesh elements [m]

% make UnitCell into Matrix
UnitRow = UnitCell(:, 3);
p = zeros(Pix); %Setting up a matrix to contain the dielectric values at positions corresponding to their mesh positions in the unit cell
for i = 1:Pix
    p(i, :) = UnitRow(((i - 1) * Pix + 1):(i * Pix)); %p equals dielectric values (p = Eps)
end
% p = rot90(p, 3);

% clear p;
% load Test2.mat

%% Setting up the concatenation of the four (Pix^2 * Pix^2) matrices
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
        % Main diagonal - Legit
        for id = 1:Pix
           Hxi{i, i}(id, id) = (1 / p(i, id)) * (2 / ax^2 + 2 / ay^2);
        end

        % Upper diagonal - Legit
        Hxi{i, i}(1, 2) = (1 / p(i, 1)^2) * (p(i, 2) - p(i, Pix) - 4 * p(i, 1)) / (4 * ay^2); 
        for id = 2:(Pix - 1)
           Hxi{i, i}(id, id + 1) = (1 / p(i, id)^2) * (p(i, id + 1) - p(i, id - 1) - 4 * p(i, id)) / (4 * ay^2); 
        end

        % Lower diagonal - Legit
        Hxi{i, i}(Pix, Pix - 1) = (1 / p(i, Pix)^2) * (-p(i, 1) + p(i, Pix - 1) - 4 * p(i, Pix)) / (4 * ay^2);
        for id = 1:(Pix - 2)
           Hxi{i, i}(id + 1, id) = (1 / p(i, id + 1)^2) * (-p(i, 2 + id) + p(i, id) - 4 * p(i, id + 1)) / (4 * ay^2);
        end
        
        % Top corner - Legit
        Hxi{i, i}(1, Pix) = exp(1i * k(m, 2) * yPeriod) * (1 / p(i, 1)^2) * (-p(i, 2) + p(i, Pix) - 4 * p(i, 1)) / (4 * ay^2);

        % Bot corner - Legit
        Hxi{i, i}(Pix, 1) = exp(-1i * k(m, 2) * yPeriod) * (1 / p(i, Pix)^2) * (p(i, 1) - p(i, Pix - 1) - 4 * p(i, Pix)) / (4 * ay^2);
     end
     
     % Upper and lower diagonal of (Pix-1) by (Pix * Pix) Matrices - Legit
     for i = 1:(Pix - 1)
         for id = 1:Pix
            Hxi{i, i + 1}(id, id) = -(1 / p(i, id)) / ax^2;
            Hxi{i + 1, i}(id, id) = -(1 / p(i + 1, id)) / ax^2;
         end
     end
     
     % Corner Pix * Pix matrix - Legit
     for id = 1:Pix
         Hxi{1, Pix}(id, id) = exp(1i * k(m, 1) * xPeriod) * (-(1 / p(1, id)) / ax^2);
         Hxi{Pix, 1}(id, id) = exp(-1i * k(m, 1) * xPeriod) * (-(1 / p(Pix, id)) / ax^2);
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

        % Main diagonal - Legit
        for id = 1:Pix
           Hyj{i, i}(id, id) = (1 / p(i, id)) * (2 / ax^2 + 2 / ay^2);
        end

        % Upper diagonal - Legit
        for id = 1:(Pix - 1)
           Hyj{i, i}(id, id + 1) = -(1 / p(i, id)) / ay^2;
        end

        % Lower diagonal - Legit
        for id = 2:Pix
           Hyj{i, i}(id, id - 1) = -(1 / p(i, id)) / ay^2;
        end
        
        % Top corner - Legit
         Hyj{i, i}(1, Pix) = exp(1i * k(m, 2) * yPeriod) * (-(1 / p(i, 1)) / ay^2);

        % Bot corner - Legit
         Hyj{i, i}(Pix, 1) = exp(-1i * k(m, 2) * yPeriod) * (-(1 / p(i, Pix)) / ay^2);
     end
     
     % Upper and lower diagonal of (Pix-1) by (Pix * Pix) Matrices - Legit
     for id = 1:Pix
         Hyj{1, 2}(id, id) = (1 / p(1, id)^2) * (p(2, id) - p(Pix, id) - 4 * p(1, id)) / (4 * ax^2);
         Hyj{Pix, Pix - 1}(id, id) = (1 / p(Pix, id)^2) * (-p(1, id) + p(Pix - 1, id) - 4 * p(Pix, id)) / (4 * ax^2);
     end
         
     for i = 2:(Pix - 1)
         for id = 1:Pix
            Hyj{i, i + 1}(id, id) = (1 / p(i, id)^2) * (p(i + 1, id) - p(i - 1, id) - 4 * p(i, id)) / (4 * ax^2); 
         end
     end
     
     for i = 1:(Pix - 2)
         for id = 1:Pix
            Hyj{i + 1, i}(id, id) = (1 / p(i + 1, id)^2) * (-p(i + 2, id) + p(i, id) - 4 * p(i + 1, id)) / (4 * ax^2);
         end
     end     
     
     % Corner Pix * Pix matrix - Legit
     for id = 1:Pix
         Hyj{1, Pix}(id, id) = exp(1i * k(m, 1) * xPeriod) * (-p(2, id) + p(Pix, id) - 4 * p(1, id)) / (4 * ax^2) * (1 / p(1, id)^2);
         Hyj{Pix, 1}(id, id) = exp(-1i * k(m, 1) * xPeriod) * (p(1, id) - p(Pix - 1, id) - 4 * p(Pix, id)) / (4 * ax^2) * (1 / p(Pix, id)^2);
     end
     
    %% Hyi (Top Right)
    % Initializes the (Pix^2 * Pix^2) matrix that contains Pix number of (Pix * Pix) matrices
    Hyi = cell(Pix);
    for i = 1:Pix
        for j = 1:Pix
            Hyi{i, j} = zeros(Pix);
        end
    end
    
    % Upper diagonal of (Pix-1) by (Pix * Pix) matrices in size - Legit
    for i = 1:(Pix - 1)
        HyiL{i, i + 1}(Pix, Pix) = p(i, 1);
        HyiR{i, i + 1}(1, 1) = p(i, Pix);
        for id = 1:(Pix - 1)
            HyiL{i, i + 1}(id, id) = p(i, id + 1);
            HyiR{i, i + 1}(id + 1, id + 1) = p(i, id);
        end
        Hyi{i, i + 1} = (-HyiL{i, i + 1} + HyiR{i, i + 1}) / (4 * ax * ay);
        for id = 1:Pix %need loop for id term in p
            Hyi{i, i + 1}(id, id) = Hyi{i, i + 1}(id, id) / p(i, id)^2; % double check the subscripts here
        end
    end
    
    % Bottom diagonal of (Pix-1) by (Pix * Pix) matrices in size - Legit
    for i = 1:(Pix - 1)
        HyiL{i + 1, i}(Pix, Pix) = p(i + 1, 1);
        HyiR{i + 1, i}(1, 1) = p(i + 1, Pix);
        for id = 1:(Pix - 1)
            HyiL{i + 1, i}(id, id) = p(i + 1, id + 1);
            HyiR{i + 1, i}(id + 1, id + 1) = p(i + 1, id);
        end
        Hyi{i + 1, i} = (HyiL{i + 1, i} - HyiR{i + 1, i}) / (4 * ax * ay);
        for id = 1:Pix
            Hyi{i + 1, i}(id, id) = Hyi{i + 1, i}(id, id) / p(i, id)^2;
        end
    end
    
    % Top right Pix * Pix matrix - Legit
    HyiL{1, Pix}(Pix, Pix) = p(1, 1);
    HyiR{1, Pix}(1, 1) = p(1, Pix);
    for id = 1:(Pix - 1)
        HyiL{1, Pix}(id, id) = p(1, id + 1);
        HyiR{1, Pix}(id + 1, id + 1) = p(1, id);
    end
    Hyi{1, Pix} = exp(1i * k(m, 1) * xPeriod) * (HyiL{1, Pix} - HyiR{1, Pix}) / (4 * ax * ay);
    for id = 1:Pix
        Hyi{1, Pix}(id, id) = Hyi{1, Pix}(id, id) / p(1, id)^2;
    end
    
    % Bottom Left Pix * Pix matrix - Legit
    HyiR{Pix, 1}(Pix, Pix) = p(Pix, 1);
    HyiL{Pix, 1}(1, 1) = p(Pix, Pix);
    for id = 1:(Pix - 1)
        HyiR{Pix, 1}(id, id) = p(Pix, id + 1);
        HyiL{Pix, 1}(id + 1, id + 1) = p(Pix, id);
    end
    Hyi{Pix, 1} = exp(-1i * k(m, 1) * xPeriod) * (-HyiL{Pix, 1} + HyiR{Pix, 1}) / (4 * ax * ay);
    for id = 1:Pix
        Hyi{Pix, 1}(id, id) = Hyi{Pix, 1}(id, id) / p(Pix, id)^2;
    end
    
    %% Hxj (Bottom Left)
    % Initializes the (Pix^2 * Pix^2) matrix that contains Pix number of (Pix * Pix) matrices
    Hxj = cell(Pix);
    for i = 1:Pix
        for j = 1:Pix
            Hxj{i, j} = zeros(Pix); % Initializes the (Pix * Pix) matrix
        end
    end
    
    % Top Right and Bottom Left corners of Pix by (Pix * Pix) matrices - Legit
    HxjR{Pix, Pix}(1, Pix) = p(1, 1);  
    HxjL{1, 1}(1, Pix) = p(Pix, 1);
    HxjL{Pix, Pix}(Pix, 1) = p(1, Pix);  
    HxjR{1, 1}(Pix, 1) = p(Pix, Pix);
    for i = 1:(Pix - 1)
        HxjR{i, i}(1, Pix) = p(i + 1, 1);
        HxjL{i + 1, i + 1}(1, Pix) = p(i, 1);    
        HxjL{i, i}(Pix, 1) = p(i + 1, Pix);
        HxjR{i + 1, i + 1}(Pix, 1) = p(i, Pix);    
    end
    for i = 1:Pix
        Hxj{i, i}(1, Pix) = (1 / p(1, Pix)^2) * exp(1i * k(m, 2) * yPeriod) * (-HxjL{i, i}(1, Pix) + HxjR{i, i}(1, Pix)) / (4 * ax * ay);
        Hxj{i, i}(Pix, 1) = (1 / p(Pix, 1)^2) * exp(-1i * k(m, 2) * yPeriod) * (-HxjL{i, i}(Pix, 1) + HxjR{i, i}(Pix, 1)) / (4 * ax * ay);
    end
    
    % Upper and bottm diagonals of Pix by (Pix * Pix) matrices - Legit
    %Bauman - None of these can be left outside of the loops because either id or i is always looping in all cases
    for i = 1:(Pix - 1) % Pix * Pix diagonal cells within the Pix^2 matrix
        for id = 1:(Pix - 1) % Selecting cell positions
            HxjL{i, i}(id, id + 1) = p(i + 1, id); % Diagonal starting at cell (1,2) in printout (left term)
            HxjR{i + 1, i + 1}(id, id + 1) = p(i, id); % Diagonal starting at cell (4,5) in printout (right term)
            HxjL{i, i}(id + 1, id) = p(i + 1, id + 1); % Diagonal starting at cell (2,1), L term
            HxjR{i + 1, i + 1}(id + 1, id) = p(i, id + 1); % Diagonal starting at cell (5,4), R term
            HxjL{Pix, Pix}(id, id + 1) = p(1, id); % Diagonal starting at cell (7,8), L term
            HxjR{1, 1}(id, id + 1) = p(Pix, id); % Diagonal starting at cell (1,2), R term
            HxjL{Pix, Pix}(id + 1, id) = p(1, id + 1); % Diagonal starting at cell (8,7), L term
            HxjR{1, 1}(id + 1, id) = p(Pix, id + 1); % Diagonal starting at cell (2,1), R term
        end
    end
    % Adds terms together - Legit
    for i = 1:Pix % For each Pix * Pix groups
        for id = 1:(Pix - 1) % Pix - 1 cells in each of these groups
            % Top diagonals:
            Hxj{i, i}(id, id + 1) = (-HxjL{i, i}(id, id + 1) + HxjR{i, i}(id, id + 1)) / (4 * ax * ay);
            % Bot diagonals:
            Hxj{i, i}(id + 1, id) = (HxjL{i, i}(id + 1, id) - HxjR{i, i}(id + 1, id)) / (4 * ax * ay);
        end
        for id = 1:Pix
            Hxj{i, i}(id, id) = Hxj{i, i}(id, id) / p(i, id)^2;
        end
    end
    
    %% Concatenating the 4 (Pix^2 * Pix^2) matrices into a 2 * (Pix^2 by Pix^2) matrix
    Hxi = cell2mat(Hxi);
    Hyi = cell2mat(Hyi);
    Hxj = cell2mat(Hxj);
    Hyj = cell2mat(Hyj);
%     H_Top{m} = cat(2, Hyj, Hyi);
%     H_Bot{m} = cat(2, Hxj, Hxi);
%     %Tests
   H_Top{m} = cat(2, Hxi, Hyi);
   H_Bot{m} = cat(2, Hxj, Hyj);
%    H_Top{m} = cat(2, Hxi, zeros(Pix^2));
%    H_Bot{m} = cat(2, zeros(Pix^2), Hyj);
   H{m} = cat(1, H_Top{m}, H_Bot{m});
%     H{m} = Hxi + Hyi + Hxj + Hyj;
% H{m} = Hxi;
end

% spy(H{1})


FDFDmatrix = H'; % Transpose (appears to do nothing to the results)
end

