function [Grid] = Quiver_complex(xfrom, xto, Nx, yfrom, yto, Ny, func)
% Setting the grid to zero value
Grid = zeros(Ny,Nx);

% Getting the x and y array
X = linspace(xfrom, xto, Nx);
Y = linspace(yfrom, yto, Ny);

% Calcuating the corresponding value for each coordinate z = x + iy
for row = 1:Nx
    for col = 1:Ny
        Grid(row,col) = func( complex( X(col), Y(row) ) );
    end
end

% Plotting the quiver
axis equal
quiver(X, Y, real(Grid), -imag(Grid),'blue');
% - imag because w = u_x - i*u_y
end