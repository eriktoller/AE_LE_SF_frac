function [Grid_1,Grid_2,Grid_t,X,Y] = Contour_major_minor(xfrom, xto, Nx, yfrom, yto, Ny, func)
% Setting the grids to zeros value
Grid_1 = zeros(Ny,Nx);
Grid_2 = zeros(Ny,Nx);
Grid_t = zeros(Ny,Nx);

% Creating the x and y array
X = linspace(xfrom, xto, Nx);
Y = linspace(yfrom, yto, Ny);

% Calculating the time for a single row
tic
for row = 1:1
    for col = 1:Ny
        [s1,s2,tp] = func( complex( X(col), Y(row) ) );
        Grid_1(row,col) = s1;
        Grid_2(row,col) = s2;
        Grid_t(row,col) = tp;
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
        [s1,s2,tp] = func( complex( X(col), Y(row) ) );
        Grid_1(row,col) = s1;
        Grid_2(row,col) = s2;
        Grid_t(row,col) = tp;
    end
    str = sprintf(formatSpec,ceil(time*(Nx-row)));
    waitbar(row/Nx,f,str)
end
close(f)

end