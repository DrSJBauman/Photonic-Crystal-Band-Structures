function [HexMatrix, radius] = CreateRegHexHole(xPix, DielectConst, rRatio)
%Creates a single unit cell for a hexagonal lattice with air holes
%surrounded by a material with the dielectric constant provided.
    
    %xPix: number of mesh elements 
    radius = rRatio * xPix;
    SF = round(2 / sqrt(3));
    
    HexMatrix = zeros(xPix, 3 * round(xPix / 2));
    Slope = sqrt(3);
    [columnsInImage, rowsInImage] = meshgrid(1:round(3 * xPix / 2), 1:xPix);
    centerX = round(3 * xPix / 4);
    centerY = round(xPix / 2);
    CirclePixels = ((rowsInImage - centerY) / (SF * radius)).^2 ...
        + ((columnsInImage - centerX) / radius).^2 <= 1;
    CirclePixels = double(CirclePixels);
    Values = cell(xPix, 1);
    Min = zeros(xPix, 1);
    Max = zeros(xPix, 1);
    for i = 1:size(HexMatrix, 1)
        Values{i} = find(CirclePixels(i, :) == 1);
        if isempty(Values{i}) == 0
            Min(i) = min(Values{i});
            Max(i) = max(Values{i});
        end
    end
    CirPixOut = zeros(xPix, 3 * round(xPix / 2));
    for i = 1:size(HexMatrix, 1)
        if Min(i) > 0
            CirPixOut(i, Min(i)) = 1;
            CirPixOut(i, Max(i)) = 1;
        end
    end
    for i = 1:size(HexMatrix, 1)
        xPoints(i) = round(i / (2 * Slope / sqrt(3)));
    end
    xPoints = xPoints';
    for i = 1:size(HexMatrix, 1)
        HexMatrix(i, xPoints(i)) = 1;
        HexMatrix(i, xPoints(i) + xPix) = 1;
    end
    HexMatrix = flipud(HexMatrix);
    HexMatrix = HexMatrix + CirPixOut;
    HexOnes = cell(xPix, 1);
    for i = 1:xPix
        HexOnes{i} = find(HexMatrix(i, :) == 1);
    end
    for i = 1:size(HexMatrix, 1)
        for j = 1:size(HexMatrix, 2)
            if j >= HexOnes{i}(1) && j <= HexOnes{i}(2)
                HexMatrix(i, j) = 1;
            end
            if j >= HexOnes{i}(size(HexOnes{i}, 2) - 1) && j <= HexOnes{i}(size(HexOnes{i}, 2))
                HexMatrix(i, j) = 1;
            end
        end
    end
    for i = 1:xPix
        M(i) = find(HexMatrix(i, :) == 1, 1);
        HexMatrix(i, :) = circshift(HexMatrix(i, :), [0 -(M(i) - 1)]);
    end
    HexMatrix = HexMatrix(1:xPix, 1:xPix);   
    for i = 1:size(HexMatrix, 1)
        for j = 1:size(HexMatrix, 2)
            if HexMatrix(i, j) == 1
               HexMatrix(i, j) = DielectConst;
            else 
               HexMatrix(i, j) = 1;
            end
        end
    end
    HexMatrix = fliplr(HexMatrix);
end

