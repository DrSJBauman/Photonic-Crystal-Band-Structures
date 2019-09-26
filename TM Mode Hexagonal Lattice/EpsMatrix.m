function [epsA] = EpsMatrix(CirclePixels)
%input: The matrix containing the circular unit cell geometry with dielectric constant values in the correct positions.
%output: A matrix containing dielectric constant values along the diagonal?

    xPix = size(CirclePixels, 1); % number of x-mesh elements (dimensionless)
    yPix = size(CirclePixels, 2); % number of y-mesh elements (dimensionless)
    
    % Initialize epsilon matrix
    epsA = zeros(xPix * yPix); % create matrix of zeros (xPix by yPix in size)
    eps = CirclePixels(1, :); % Vector eps equals CirclePixels 1st row all columns
    
    % Create epsilon vector
    for i = 2:yPix
        eps = [eps, CirclePixels(i, :)]; %Is eps a vector or a matrix?
    end
    
    % Create epsilon matrix
    for i = 1:(xPix * yPix)
        epsA(i, i) = eps(i); %Every position along the diagonal of epsA is assigned value from eps vector
    end
%epsA = flipud(epsA);
%epsA = fliplr(epsA);
end

