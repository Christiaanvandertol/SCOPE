function int = Sint(y,x)

    % Simpson integration
    % x and y must be any vectors (rows, columns), but of the same length
    % x must be a monotonically increasing series
    
    % WV Jan. 2013, for SCOPE 1.40
    
    nx   = length(x);
    if size(x,1) == 1
        x = x';
    end
    if size(y,1) ~= 1
        y = y';
    end
    step = x(2:nx) - x(1:nx-1);
    mean = .5 * (y(1:nx-1) + y(2:nx));
    int  = mean * step;
end