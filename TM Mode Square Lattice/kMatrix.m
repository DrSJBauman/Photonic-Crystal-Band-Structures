function [FDFDmatrix] = kMatrix(xPeriod, yPeriod, xPix, yPix, k)
    % Initialize tauOneDiagUnitCell cell array
    tauOneDiagUnitCell = cell(size(k, 1), 1);
    % Matrix denominators
    Ax = xPeriod / xPix;
    Ay = yPeriod / yPix;
    % Initialize tauOneDiagUnitCell matrix
    for m = 1:size(k, 1)
        tauOneDiagUnitCell{m} = zeros(xPix);
    end
    % Create diagonal along tauOneDiagUnitCell matrix
    for m = 1:size(k, 1)
        for j = 1:xPix
            tauOneDiagUnitCell{m}(j, j) = (2 / Ax^2) + (2 / Ay^2);
        end
    end
    % Create adjacent diagonals
    for m = 1:size(k, 1)
        for j = 1:(xPix - 1)
            tauOneDiagUnitCell{m}(j, (j + 1)) =  - (1 / Ay^2);
        end
    end
    for m = 1:size(k, 1)
        for j = 2:xPix
            tauOneDiagUnitCell{m}(j, (j - 1)) =  - (1 / Ay^2);
        end
    end
    % Create outer points
    for m = 1:size(k, 1)
        tauOneDiagUnitCell{m}(1, xPix) = exp(1i * k(m, 2) * yPeriod) * ( - (1 / Ay^2));
        tauOneDiagUnitCell{m}(xPix, 1) = exp(-1i * k(m, 2) * yPeriod) * ( - (1 / Ay^2));
    end

    % Initialize tauTwoDiagUnitCell and tauTwoDiagUnitCellConj
    tauTwoDiagUnitCell = cell(size(k, 1), 1);
    tauTwoDiagUnitCellConj = cell(size(k, 1), 1);
    % Create tauTwoDiagUnitCell and tauTwoDiagUnitCellConj
    for m = 1:size(k, 1)
        for j = 1:xPix
            tauTwoDiagUnitCell{m}(j, j) =   - (1 / Ax^2);
            tauTwoDiagUnitCellConj{m}(j, j) = conj(tauTwoDiagUnitCell{m}(j, j));
        end
    end

    % Initialize tauThreeDiagUnitCell and tauThreeDiagUnitCellConj
    tauThreeDiagUnitCell = cell(size(k, 1), 1);
    tauThreeDiagUnitCellConj = cell(size(k, 1), 1);
    % Create tauThreeDiagUnitCell and tauThreeDiagUnitCellConj
    for m = 1:size(k, 1)
        for j = 1:xPix
            tauThreeDiagUnitCell{m}(j, j) = exp(1i * k(m, 1) * yPeriod) * ( - (1 / Ax^2));
            tauThreeDiagUnitCellConj{m}(j, j) = conj(tauThreeDiagUnitCell{m}(j, j));
        end
    end

    % Initialize array of final matrix
    FDFDmatrix = cell(size(k, 1), 1);
    % Create empty unit cell matrix
    zeroUnitCell = zeros(xPix);
    % Initialize final matrix
    for m = 1:size(k, 1)
        for j = 1:yPix
            for g = 1:xPix
                FDFDmatrix{m}{j, g} = zeroUnitCell;
            end
        end
    end
    % Create cell diagonals
    for m = 1:size(k, 1)
        for j = 1:yPix
            FDFDmatrix{m}{j, j} = tauOneDiagUnitCell{m};
        end
    end
    % Create adjacent cell diagonals
    for m = 1:size(k, 1)
        for j = 1:(yPix - 1)
            FDFDmatrix{m}{j, (j + 1)} = tauTwoDiagUnitCell{m};
        end
    end
    for m = 1:size(k, 1)
        for j = 2:yPix
            FDFDmatrix{m}{j, (j - 1)} = tauTwoDiagUnitCellConj{m};
        end
    end
    % Create outer matrix
    for m = 1:size(k, 1)
        FDFDmatrix{m}{1, yPix} = tauThreeDiagUnitCell{m};
        FDFDmatrix{m}{yPix, 1} = tauThreeDiagUnitCellConj{m};
    end
    % Convert final matrix cell array to array
    for m = 1:size(k, 1)
        FDFDmatrix{m} = cell2mat(FDFDmatrix{m});
    end
end

