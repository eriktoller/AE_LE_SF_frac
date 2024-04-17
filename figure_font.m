function [] = figure_font(font_size)
%FIGURE_FONT Set the font properties for figures
%   The function sets the font, font size and interpreter for the
%   figures (axes, legends and title).
%
%   VARIBALES
%   font_size - font size for text in figure
%
%   LATEST UPDATE
%   2023-01-15
%
%   AUTHOR
%   Erik Toller,
%   Department of Earth Sciences, Uppsala University, Sweden

switch nargin
    case 1 % if font size is given
        set(groot,'defaultAxesTickLabelInterpreter','latex');
        set(groot,'defaulttextinterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');
        set(groot,'defaultAxesFontSize',font_size)
    otherwise % if no font size is given
        set(groot,'defaultAxesTickLabelInterpreter','latex');
        set(groot,'defaulttextinterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');
        set(groot,'defaultAxesFontSize',15)
end

end

