function [CirclePixels, radius] = CreateHole(xPix, yPix, DConst, rRatio)
%Creates a single unit cell for a square lattice with air holes
%surrounded by a material with the dielectric constant provided.
    radius = rRatio * xPix;
%     radius = round((xPix / xPeriod) * radius);
    [columnsInImage, rowsInImage] = meshgrid(1:xPix, 1:yPix);
    centerX = xPix / 2;
    centerY = yPix / 2;
    CirclePixels = (rowsInImage - centerY).^2 ...
        + (columnsInImage - centerX).^2 <= radius.^2;
    CirclePixels = double(CirclePixels);
    for i = 1:size(CirclePixels, 1)
        for j = 1:size(CirclePixels, 2)
            if CirclePixels(i, j) == 0
                CirclePixels(i, j) = DConst; %Dielectric constant for surrounding material
            else 
                CirclePixels(i, j) = 1;        %Dielectric constant of the holes 
            end
        end
    end
end

