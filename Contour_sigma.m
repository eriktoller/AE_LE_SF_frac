function [Grid_11,Grid_22,Grid_12,X,Y] = Contour_sigma(xfrom, xto, Nx, yfrom, yto, Ny, func)
% Setting the grids to zeros value
Grid_11 = zeros(Ny,Nx);
Grid_22 = zeros(Ny,Nx);
Grid_12 = zeros(Ny,Nx);

% Creating the x and y array
X = linspace(xfrom, xto, Nx);
Y = linspace(yfrom, yto, Ny);

% Calculating the time for a single row
tic
for row = 1:1
    for col = 1:Ny
        [s1,s2,s3] = func( complex( X(col), Y(row) ) );
        Grid_11(row,col) = s1;
        Grid_22(row,col) = s2;
        Grid_12(row,col) = s3;
    end
end
% Creating the progress bar and estimating the time left
time = toc;
formatSpec = 'Time remining: %d sec';
str = sprintf(formatSpec,ceil(time*(Nx-row)));
f = waitbar(0,str,'Name','Plotting contour');

% Calcuating for the remaining rows
for row = 2:Nx
    for col = 1:Ny
        [s1,s2,s3] = func( complex( X(col), Y(row) ) );
        Grid_11(row,col) = s1;
        Grid_22(row,col) = s2;
        Grid_12(row,col) = s3;
    end
    str = sprintf(formatSpec,ceil(time*(Nx-row)));
    waitbar(row/Nx,f,str)
end
close(f)

end