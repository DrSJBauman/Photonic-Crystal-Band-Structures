function [CirclePixels, radius] = CreatePost(xPix, yPix, DConst, rRatio)
%Creates a single unit cell for a square lattice of posts with the
%dielectric constant provided, surrounded by air.
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
                CirclePixels(i, j) = 1; %Dielectric constant for surrounding material
            else 
                CirclePixels(i, j) = DConst; %Dielectric constant of the posts 
            end
        end
    end
end

