function [TEtrix] = TEMatrix(xPeriod, yPeriod, xPix, yPix, k, UnitCell)
    %% Parameters:
    Pix = xPix;
    UnitRow = 1 ./ UnitCell(:, 3);
    f = zeros(Pix); %Setting up a matrix to contain the dielectric values at positions corresponding to their mesh positions in the unit cell
    for i = 1:Pix
        f(i, :) = UnitRow(((i - 1) * Pix + 1):(i * Pix));
    end
    Nmesh = xPix * yPix;
    %% Create TL Matrix
    % Initialize matrix cell array
    TEtrixDiag = cell(size(k, 1), xPix);
    % Matrix denominators
    ALx = xPeriod / xPix;
    ALy = yPeriod / yPix;
    % Initialize tauOneDiagUnitCell matrix
    for m = 1:size(k, 1)
        for i = 1:xPix
            TEtrixDiag{m, i} = zeros(xPix);
        end
    end
    % Creates diagonal elements of tauOneDiagUnitCell matrix
    for m = 1:size(k, 1)
        for i = 1:xPix
            for j = 1:xPix
                TEtrixDiag{m, i}(j, j) = f(i, j) * ((2 / ALx^2) + (2 / ALy^2));
            end
        end
    end
    % Creates left diagonal elements of tauOneDiagUnitCell matrix
    for m = 1:size(k, 1)
        for i = 1:xPix
            for j = 2:xPix
                if j + 1 == xPix + 1
                    TEtrixDiag{m, i}(j, j - 1) = (1 / (4 * ALy^2)) * (f(i, 1) - f(i, j - 1) - 4 * f(i, j)); 
                else
                    TEtrixDiag{m, i}(j, j - 1) = (1 / (4 * ALy^2)) * (f(i, j + 1) - f(i, j - 1) - 4 * f(i, j));
                end
            end
        end
    end
    % Creates right diagonal elements of tauOneDiagUnitCell matrix
    for m = 1:size(k, 1)
        for i = 1:xPix
            for j = 2:xPix
                if j + 1 == xPix + 1
                    TEtrixDiag{m, i}(j - 1, j) = (1 / (4 * ALy^2)) * (-f(i, j) + f(i, 1) - 4 * f(i, j - 1)); 
                else
                    TEtrixDiag{m, i}(j - 1, j) = (1 / (4 * ALy^2)) * (-f(i, j) + f(i, j + 1) - 4 * f(i, j - 1));
                end
            end
        end
    end
    % Creates corner elements of tauOneDiagUnitCell matrix
    for m = 1:size(k, 1)
        for j = 1:xPix
            TEtrixDiag{m, j}(1, xPix) = (exp(1i * k(m, 2) * yPeriod) / (4 * ALy^2)) * (f(j, 2) - f(j, xPix) - 4 * f(j, 1));
            TEtrixDiag{m, j}(xPix, 1) = (exp(-1i * k(m, 2) * yPeriod) / (4 * ALy^2)) * (-f(j, 1) + f(j, (xPix - 1)) - 4 * f(j, xPix));
        end
    end

    % Initialize TEtrixDiagTop matrix
    TEtrixDiagTop = cell(size(k, 1), (xPix - 1));
    for m = 1:size(k, 1)
        for i = 1:(xPix - 1)
            TEtrixDiagTop{m, i} = zeros(xPix);
        end
    end
    % Create TEtrixDiagTop matrix
    for m = 1:size(k, 1)
        for i = 1:(xPix - 1)
            for j = 1:xPix
                TEtrixDiagTop{m, i}(j, j) = -f(i, j) / ALx^2;
            end
        end
    end

    % Initialize TEtrixDiagBot matrix
    TEtrixDiagBot = cell(size(k, 1), (xPix - 1));
    for m = 1:size(k, 1)
        for i = 1:(xPix - 1)
            TEtrixDiagBot{m, i} = zeros(xPix);
        end
    end
    % Create TEtrixDiagBot matrix
    for m = 1:size(k, 1)
        for i = 1:(xPix - 1)
            for j = 1:xPix     
                TEtrixDiagBot{m, i}(j, j) = -f(i + 1, j) / ALx^2;
            end
        end
    end

    % Initialize TEtrixDiagRight matrix
    TEtrixDiagRight = cell(size(k, 1), 1);
    for m = 1:size(k, 1)
        TEtrixDiagRight{m} = zeros(xPix);
    end
    % Create TEtrixDiagRight matrix
    for m = 1:size(k, 1)
        for j = 1:xPix
            TEtrixDiagRight{m}(j, j) = -(exp(1i * k(m, 1) * xPeriod) / ALx^2) * f(1, j);
        end
    end

    % Initialize TEtrixDiagLeft matrix
    TEtrixDiagLeft = cell(size(k, 1), 1);
    for m = 1:size(k, 1)
        TEtrixDiagLeft{m} = zeros(xPix);
    end
    % Create TEtrixDiagLeft matrix
    for m = 1:size(k, 1)
        for j = 1:xPix
            TEtrixDiagLeft{m}(j, j) = -(exp(-1i * k(m, 1) * xPeriod) / ALx^2) * f(xPix, j);
        end
    end

    % Initialize Final cell Matrix
    TEtrixCell = cell(xPix);
    TEtrixZeros = zeros(xPix);
    for m = 1:size(k, 1)
        for i = 1:xPix
            for j = 1:xPix
                TEtrixCell{m}{i, j} = TEtrixZeros;
            end
        end
    end
    % Create diagonal entrees
    for m = 1:size(k, 1)
        for i = 1:xPix
            TEtrixCell{m}{i, i} = TEtrixDiag{m, i};
        end
    end
    % Create top diagonal entrees
    for m = 1:size(k, 1)
        for i = 2:xPix
            TEtrixCell{m}{i - 1, i} = TEtrixDiagTop{m, i - 1};
        end
    end
    % Create bot diagonal entrees
    for m = 1:size(k, 1)
        for i = 2:xPix
            TEtrixCell{m}{i, i - 1} = TEtrixDiagBot{m, i - 1};
        end
    end
    % Create corner entrees
    for m = 1:size(k, 1)
        TEtrixCell{m}{1, xPix} = TEtrixDiagRight{m};
        TEtrixCell{m}{xPix, 1} = TEtrixDiagLeft{m};
    end
    % Initialize final matrix
    TEtrixTL = cell(size(k, 1), 1);
    for m = 1:size(k, 1)    
        TEtrixTL{m} = zeros(Nmesh);
    end
    % Convert cell array to double array
    for m = 1:size(k, 1)
        TEtrixTL{m} = cell2mat(TEtrixCell{m});
    end

    %% Create BR Matrix
    % Initialize matrix cell array
    TEtrixDiagBR = cell(size(k, 1), xPix);
    % Initialize tauOneDiagUnitCell matrix
    for m = 1:size(k, 1)
        for i = 1:xPix
            TEtrixDiagBR{m, i} = zeros(xPix);
        end
    end
    % Creates diagonal elements of tauOneDiagUnitCell matrix
    for m = 1:size(k, 1)
        for i = 1:xPix
            for j = 1:xPix
                TEtrixDiagBR{m, i}(j, j) = f(i, j) * ((2 / ALx^2) + (2 / ALy^2));
            end
        end
    end
    % Creates left diagonal elements of tauOneDiagUnitCell matrix
    for m = 1:size(k, 1)
        for i = 1:xPix
            for j = 2:xPix
                TEtrixDiagBR{m, i}(j, j - 1) = -f(i, j) / ALy^2; 
            end
        end
    end
    % Creates right diagonal elements of tauOneDiagUnitCell matrix
    for m = 1:size(k, 1)
        for i = 1:xPix
            for j = 2:xPix
                TEtrixDiagBR{m, i}(j - 1, j) = -f(i, j - 1) / ALy^2; 
            end
        end
    end
    % Creates corner elements of tauOneDiagUnitCell matrix
    for m = 1:size(k, 1)
        for j = 1:xPix
            TEtrixDiagBR{m, j}(1, xPix) = -(exp(1i * k(m, 2) * yPeriod) / ALy^2) * f(i, 1);
            TEtrixDiagBR{m, j}(xPix, 1) = -(exp(-1i * k(m, 2) * yPeriod) / ALy^2) * f(i, xPix);
        end
    end

    % Initialize TEtrixDiagTop matrix
    TEtrixDiagTopBR = cell(size(k, 1), (xPix - 1));
    for m = 1:size(k, 1)
        for i = 1:(xPix - 1)
            TEtrixDiagTopBR{m, i} = zeros(xPix);
        end
    end
    % Create TEtrixDiagTop matrix
    for m = 1:size(k, 1)
        for i = 2:(xPix - 1)
            for j = 1:xPix
                TEtrixDiagTopBR{m, 1}(j, j) = (1 / 4 * ALx^2) * (-f(2, j) + f(xPix, j) - 4 * f(1, j));
                TEtrixDiagTopBR{m, i}(j, j) = (1 / 4 * ALx^2) * (-f((i + 1), j) + f(i - 1, j) - 4 * f(i, j));
            end
        end
    end

    % Initialize TEtrixDiagBot matrix
    TEtrixDiagBotBR = cell(size(k, 1), (xPix - 1));
    for m = 1:size(k, 1)
        for i = 1:(xPix - 1)
            TEtrixDiagBotBR{m, i} = zeros(xPix);
        end
    end
    % Create TEtrixDiagBot matrix
    for m = 1:size(k, 1)
        for i = 2:(xPix - 1)
            for j = 1:xPix     
                TEtrixDiagBotBR{m, 1}(j, j) = (1 / 4 * ALx^2) * (f(xPix, j) - f(1, j) - 4 * f(2, j));
                TEtrixDiagBotBR{m, i}(j, j) = (1 / 4 * ALx^2) * (f(i - 1, j) - f(i, j) - 4 * f(i + 1, j));
            end
        end
    end

    % Initialize TEtrixDiagRight matrix
    TEtrixDiagRightBR = cell(size(k, 1), 1);
    for m = 1:size(k, 1)
        TEtrixDiagRightBR{m} = zeros(xPix);
    end
    % Create TEtrixDiagRight matrix
    for m = 1:size(k, 1)
        for j = 1:xPix
            TEtrixDiagRightBR{m}(j, j) = (exp(1i * k(m, 1) * xPeriod) / (4 * ALx^2)) * (f(2, j) - f(xPix, j) - 4 * f(1, j));
        end
    end

    % Initialize TEtrixDiagLeft matrix
    TEtrixDiagLeftBR = cell(size(k, 1), 1);
    for m = 1:size(k, 1)
        TEtrixDiagLeftBR{m} = zeros(xPix);
    end
    % Create TEtrixDiagLeft matrix
    for m = 1:size(k, 1)
        for j = 1:xPix
            TEtrixDiagLeftBR{m}(j, j) = (exp(-1i * k(m, 1) * xPeriod) / (4 * ALx^2)) * (-f(1, j) + f(2, j) - 4 * f(xPix, j));
        end
    end

    % Initialize Final cell Matrix
    TEtrixCellBR = cell(xPix);
    for m = 1:size(k, 1)
        for i = 1:xPix
            for j = 1:xPix
                TEtrixCellBR{m}{i, j} = TEtrixZeros;
            end
        end
    end
    % Create diagonal entrees
    for m = 1:size(k, 1)
        for i = 1:xPix
            TEtrixCellBR{m}{i, i} = TEtrixDiagBR{m, i};
        end
    end
    % Create top diagonal entrees
    for m = 1:size(k, 1)
        for i = 2:xPix
            TEtrixCellBR{m}{i - 1, i} = TEtrixDiagTopBR{m, i - 1};
        end
    end
    % Create bot diagonal entrees
    for m = 1:size(k, 1)
        for i = 2:xPix
            TEtrixCellBR{m}{i, i - 1} = TEtrixDiagBotBR{m, i - 1};
        end
    end
    % Create corner entrees
    for m = 1:size(k, 1)
        TEtrixCellBR{m}{1, xPix} = TEtrixDiagRightBR{m};
        TEtrixCellBR{m}{xPix, 1} = TEtrixDiagLeftBR{m};
    end
    % Initialize final matrix
    TEtrixBR = cell(size(k, 1), 1);
    for m = 1:size(k, 1)    
        TEtrixBR{m} = zeros(Nmesh);
    end
    % Convert cell array to double array
    for m = 1:size(k, 1)
        TEtrixBR{m} = cell2mat(TEtrixCellBR{m});
    end
    
    %% Create TR Matrix
    % Initialize matrix cell array
    TEtrixDiagTopTR = cell(size(k, 1), xPix - 1);
    % Initialize DiagTop matrix
    for m = 1:size(k, 1)
        for i = 1:(xPix - 1)
            TEtrixDiagTopTR{m, i} = zeros(xPix);
        end
    end
    % Create DiagTop matrix
    for m = 1:size(k, 1)
        for i = 1:(xPix - 1)
            TEtrixDiagTopTR{m, i}(1, 1) = (1 / (4 * ALx * ALy)) * (f(i, 2) - f(i, xPix));
            TEtrixDiagTopTR{m, i}(2, 2) = (1 / (4 * ALx * ALy)) * (f(i, 3) - f(i, 1));
            for j = 3:xPix
            TEtrixDiagTopTR{m, i}(j, j) = (1 / (4 * ALx * ALy)) * (f(i, j - 2) - f(i, j - 1)); 
            end
        end
    end
    % Initialize matrix cell array
    TEtrixDiagRightTR = cell(size(k, 1), 1);
    
    % Initialize DiagRight matrix
    for m = 1:size(k, 1)
        TEtrixDiagRightTR{m} = zeros(xPix);
    end   
    
    % Create DiagRight matrix
    for m = 1:size(k, 1)
        TEtrixDiagRightTR{m}(1, 1) = (exp(1i * k(m, 1) * xPeriod) / (4 * ALx * ALy)) * (-f(1, 2) + f(1, xPix));
        TEtrixDiagRightTR{m}(2, 2) = (exp(1i * k(m, 1) * xPeriod) / (4 * ALx * ALy)) * (-f(1, 3) + f(1, 1));
        for j = 3:xPix
            TEtrixDiagRightTR{m}(j, j) = (exp(1i * k(m, 1) * xPeriod) / (4 * ALx * ALy)) * (-f(1, j - 2) + f(1, j - 1));
        end
    end
    
    % Initialize matrix cell array
    TEtrixDiagBotTR = cell(size(k, 1), xPix - 1);
    % Initialize DiagBot matrix
    for m = 1:size(k, 1)
        for i = 1:(xPix - 1)
            TEtrixDiagBotTR{m, i} = zeros(xPix);
        end
    end
    
    % Create DiagBot matrix
    for m = 1:size(k, 1)
        for i = 1:(xPix - 1)
            TEtrixDiagBotTR{m, i}(1, 1) = (1 / (4 * ALx * ALy)) * (-f(i + 1, 2) + f(i + 1, xPix));
            TEtrixDiagBotTR{m, i}(2, 2) = (1 / (4 * ALx * ALy)) * (-f(i + 1, 3) + f(i + 1, 1));
            for j = 3:xPix
                TEtrixDiagBotTR{m, i}(j, j) = (1 / (4 * ALx * ALy)) * (-f(i + 1, j - 2) + f(i + 1, j - 1)); 
            end
        end
    end
    
    % Initialize matrix cell array
    TEtrixDiagLeftTR = cell(size(k, 1), 1);
    
   % Initialize DiagLeft matrix
    for m = 1:size(k, 1)
        TEtrixDiagLeftTR{m} = zeros(xPix);
    end   
    
    % Create DiagLeft matrix
    for m = 1:size(k, 1)
        TEtrixDiagLeftTR{m}(1, 1) = (exp(-1i * k(m, 1) * xPeriod) / (4 * ALx * ALy)) * (f(xPix, 2) - f(xPix, xPix));
        TEtrixDiagLeftTR{m}(2, 2) = (exp(-1i * k(m, 1) * xPeriod) / (4 * ALx * ALy)) * (f(xPix, 3) - f(xPix, 1));
        for j = 3:xPix
            TEtrixDiagLeftTR{m}(j, j) = (exp(-1i * k(m, 1) * xPeriod) / (4 * ALx * ALy)) * (f(xPix, j - 2) - f(xPix, j - 1));
        end
    end  
    
    % Initialize Final cell Matrix
    TEtrixCellTR = cell(xPix);
    for m = 1:size(k, 1)
        for i = 1:xPix
            for j = 1:xPix
                TEtrixCellTR{m}{i, j} = TEtrixZeros;
            end
        end
    end
    
    % Create Final TEtrixCellTR Matrix
    for m = 1:size(k, 1)
        TEtrixCellTR{m}{1, xPix} = TEtrixDiagRightTR{m};
        TEtrixCellTR{m}{xPix, 1} = TEtrixDiagLeftTR{m};
        for j = 2:xPix
            TEtrixCellTR{m}{j - 1, j} = TEtrixDiagTopTR{m, j - 1};
            TEtrixCellTR{m}{j, j - 1} = TEtrixDiagBotTR{m, j - 1};
        end
    end
    
    % Initialize final matrix
    TEtrixTR = cell(size(k, 1), 1);
    for m = 1:size(k, 1)    
        TEtrixTR{m} = zeros(Nmesh);
    end
    % Convert cell array to double array
    for m = 1:size(k, 1)
        TEtrixTR{m} = cell2mat(TEtrixCellTR{m});
    end
    
    %% Create BL Matrix
    % Initialize matrix cell array
    TauBLR = cell(size(k, 1), xPix);
    % Initialize TauBL matrix
    for m = 1:size(k, 1)
        for i = 1:xPix
            TauBLR{m, i} = zeros(xPix);
        end
    end
    % Create TauBL matrix
    for m = 1:size(k, 1)
        for i = 3:xPix
            for j = 2:xPix
                TauBLR{m, 1}(j - 1, j) = (1 / (4 * ALx * ALy)) * (-f(xPix, j - 1));
                TauBLR{m, 2}(j - 1, j) = (1 / (4 * ALx * ALy)) * (-f(1, j - 1));
                TauBLR{m, i}(j - 1, j) = (1 / (4 * ALx * ALy)) * (-f(i - 1, j - 1));
            end
        end
    end
    for m = 1:size(k, 1)
        for i = 3:xPix
            for j = 2:xPix
                TauBLR{m, 1}(j, j - 1) = (1 / (4 * ALx * ALy)) * (f(xPix, j));
                TauBLR{m, 2}(j, j - 1) = (1 / (4 * ALx * ALy)) * (f(1, j));
                TauBLR{m, i}(j, j - 1) = (1 / (4 * ALx * ALy)) * (f(i - 1, j));
            end
        end
    end

    % Initialize matrix cell array
    TauBLL = cell(size(k, 1), xPix);
    % Initialize TauBL matrix
    for m = 1:size(k, 1)
        for i = 1:xPix
            TauBLL{m, i} = zeros(xPix);
        end
    end
    % Create TauBL matrix
    for m = 1:size(k, 1)
        for i = 1:(xPix - 1)
            for j = 2:xPix
                TauBLL{m, i}(j - 1, j) = (1 / (4 * ALx * ALy)) * (f(i + 1, j - 1));
                TauBLL{m, xPix}(j - 1, j) = (1 / (4 * ALx * ALy)) * (f(1, j - 1));
            end
        end
    end
    for m = 1:size(k, 1)
        for i = 1:(xPix - 1)
            for j = 2:xPix
                TauBLL{m, i}(j, j - 1) = (1 / (4 * ALx * ALy)) * (-f(i + 1, j));
                TauBLL{m, xPix}(j, j - 1) = (1 / (4 * ALx * ALy)) * (-f(1, j));
            end
        end
    end
    
    % Initialize matrix cell array
    TauBL = cell(size(k, 1), xPix);
    % Initialize TauBL matrix
    for m = 1:size(k, 1)
        for i = 1:xPix
            TauBL{m, i} = TauBLR{m, i} + TauBLL{m, i};
        end
    end
    
    TauCR = cell(size(k, 1), xPix);
    for m = 1:size(k, 1)
        for i = 1:xPix
            TauCR{m, i} = zeros(xPix);
        end
    end
    for m = 1:size(k, 1)
        for i = 3:xPix
            TauCR{m, 1}(1, xPix) = (exp(1i * k(m, 2) * yPeriod) / (4 * ALx * ALy)) * f(xPix, 1);
            TauCR{m, 2}(1, xPix) = (exp(1i * k(m, 2) * yPeriod) / (4 * ALx * ALy)) * f(1, 1);
            TauCR{m, i}(1, xPix) = (exp(1i * k(m, 2) * yPeriod) / (4 * ALx * ALy)) * f(i - 1, 1);
        end
    end
                
    TauCL = cell(size(k, 1), xPix);
    for m = 1:size(k, 1)
        for i = 1:xPix
            TauCL{m, i} = zeros(xPix);
        end
    end
    for m = 1:size(k, 1)
        for i = 1:(xPix - 1)
            TauCL{m, i}(1, xPix) = (exp(1i * k(m, 2) * yPeriod) / (4 * ALx * ALy)) * (-f(i + 1, 1));
            TauCL{m, xPix}(1, xPix) = (exp(1i * k(m, 2) * yPeriod) / (4 * ALx * ALy)) * (-f(1, 1));
        end
    end
    TauC = cell(size(k, 1), xPix);
    for m = 1:size(k, 1)
        for i = 1:xPix
            TauC{m, i} = TauCR{m, i} + TauCL{m, i};
        end
    end
    
    TauBCR = cell(size(k, 1), xPix);
    for m = 1:size(k, 1)
        for i = 1:xPix
            TauBCR{m, i} = zeros(xPix);
        end
    end
    for m = 1:size(k, 1)
        for i = 3:xPix
            TauBCR{m, 1}(1, xPix) = (exp(-1i * k(m, 2) * yPeriod) / (4 * ALx * ALy)) * (-f(xPix, xPix));
            TauBCR{m, 2}(1, xPix) = (exp(-1i * k(m, 2) * yPeriod) / (4 * ALx * ALy)) * (-f(1, xPix));
            TauBCR{m, i}(1, xPix) = (exp(-1i * k(m, 2) * yPeriod) / (4 * ALx * ALy)) * (-f(i - 1, xPix));
        end
    end
                
    TauBCL = cell(size(k, 1), xPix);
    for m = 1:size(k, 1)
        for i = 1:xPix
            TauBCL{m, i} = zeros(xPix);
        end
    end
    for m = 1:size(k, 1)
        for i = 1:(xPix - 1)
            TauBCL{m, i}(1, xPix) = (exp(-1i * k(m, 2) * yPeriod) / (4 * ALx * ALy)) * f(i + 1, xPix);
            TauBCL{m, xPix}(1, xPix) = (exp(-1i * k(m, 2) * yPeriod) / (4 * ALx * ALy)) * f(1, xPix);
        end
    end
    TauBC = cell(size(k, 1), xPix);
    for m = 1:size(k, 1)
        for i = 1:xPix
            TauBC{m, i} = TauBCR{m, i} + TauBCL{m, i};
        end
    end
    
    Tau = cell(size(k, 1), xPix);
    for m = 1:size(k, 1)
        for i = 1:xPix
            Tau{m, i} = TauBL{m, i} + TauC{m, i} + TauBC{m, i};
        end
    end
    
    % Initialize Final cell Matrix
    TEtrixCellBL = cell(xPix);
    for m = 1:size(k, 1)
        for i = 1:xPix
            for j = 1:xPix
                TEtrixCellBL{m}{i, j} = TEtrixZeros;
            end
        end
    end
    
    % Create Final TEtrixCellTR Matrix
    for m = 1:size(k, 1)
        for i = 1:xPix
            TEtrixCellBL{m}{i, i} = Tau{m, i};
        end
    end
    
    % Initialize final matrix
    TEtrixBL = cell(size(k, 1), 1);
    for m = 1:size(k, 1)    
        TEtrixBL{m} = zeros(Nmesh);
    end
    % Convert cell array to double array
    for m = 1:size(k, 1)
        TEtrixBL{m} = cell2mat(TEtrixCellBL{m});
    end
    
    TEtrixTop = cell(size(k, 1), 1);
    TEtrixBot = cell(size(k, 1), 1);
    TEtrix = cell(size(k, 1), 1);
    for m = 1:size(k, 1)
        TEtrixTop{m} = cat(2, TEtrixTL{m}, TEtrixTR{m});
        TEtrixBot{m} = cat(2, TEtrixBL{m}, TEtrixBR{m});
        TEtrix{m} = cat(1, TEtrixTop{m}, TEtrixBot{m});
    end
end

