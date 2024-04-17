function [Grid,X,Y] = Contour_flownet(xfrom, xto, Nx, yfrom, yto, Ny, func)
% Set the grid to matrix of zero
Grid = zeros(Ny,Nx);
% Get the X and Y vectors
X = linspace(xfrom, xto, Nx);
Y = linspace(yfrom, yto, Ny);

% Clock the first iteration
tic
for row = 1:1
    for col = 1:Ny
        Grid(row,col) = func( complex( X(col), Y(row) ) );
    end
end
time = toc;
% Create the progress bar and estimate the time
formatSpec = 'Time remining: %d sec';
str = sprintf(formatSpec,ceil(time*(Nx-row)));
f = waitbar(0,str,'Name','Plotting contour');

% Calcuate the remaining rows and update the progress bar
for row = 2:Nx
    for col = 1:Ny
        Grid(row,col) = func( complex( X(col), Y(row) ) );
    end
    str = sprintf(formatSpec,ceil(time*(Nx-row)));
    waitbar(row/Nx,f,str)
end

% Close the progress bar
close(f)
end