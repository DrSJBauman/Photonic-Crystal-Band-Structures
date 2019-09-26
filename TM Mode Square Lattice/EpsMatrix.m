function [epsA] = EpsMatrix(CirclePixels)
    xPix = size(CirclePixels, 1); % number of x-mesh elements (dimensionless)
    yPix = size(CirclePixels, 2);
    % Initialize epsilon matrix
    epsA = zeros(xPix * yPix);
    eps = CirclePixels(1, :);
    % Create epsilon vector
    for i = 2:yPix
        eps = [eps, CirclePixels(i, :)];
    end
    % Create epsilon matrix
    for i = 1:(xPix * yPix)
        epsA(i, i) = eps(i);
    end
end

