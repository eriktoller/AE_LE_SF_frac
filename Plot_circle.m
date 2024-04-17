function [] = Plot_circle(z,r,plot_prop,ln_wdth)
% PLOT_CIRCLE Program which plots a line between complex points
%   The program plots a line between the complex points given in a vector.
%
%   VARIABLES
%   z - center point of circles (complex vector)
%   r - radius of circles (double vector)
%   plot_prop - plotting properties, e.g. 'black .-' (string)
%   ln_wdth - line width (double)
%
%   LATEST UPDATE
%   2021-10-26
%
%   AUTHOR
%   Erik Toller,
%   Department of Earth Sciences, Uppsala University, Sweden

% Plot circle with given radius at each given center point
for ii = 1:length(z)
    % Create 100 evenly distributed points along the circle
    theta = linspace(0,2*pi,100);
    z_temp = exp(1i.*theta).*r(ii) + z(ii);
    
    % Create x and y vecotrs from the complex points
    x = real(z_temp);
    y = imag(z_temp);
    
    % Plot the circle
    plot(x,y,plot_prop,'LineWidth',ln_wdth)
end
end

