function [] = Plot_line(zfrom,zto,plot_prop,ln_wdth)
% PLOT_LINE Program which plots a line between complex points
%   The program plots a line between the complex points given in a vector.
%
%   VARIABLES
%   zfrom - start point of line (complex vector)
%   zto - end point of line (complex vector)
%   plot_prop - plotting properties, e.g. 'black .-' (string)
%   ln_wdth - line width (double)
%
%   LATEST UPDATE
%   2021-10-26
%
%   AUTHOR
%   Erik Toller,
%   Department of Earth Sciences, Uppsala University, Sweden

% Plot one line between each set of points in the vectors
for ii = 1:length(zfrom)
    % Plotting a line between two points
    x = [real(zfrom(ii)), real(zto(ii))];
    y = [imag(zfrom(ii)), imag(zto(ii))];
    plot(x,y,plot_prop,'LineWidth',ln_wdth)
end
end
