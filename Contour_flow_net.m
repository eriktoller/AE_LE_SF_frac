function [M1,M2,lvs_im,lvs_re] = Contour_flow_net(X,Y,Grid,lvs,ln_wdth)
% CONTOUR_FLOW_NET Program which contours the flow net for a given grid
%   The program contours a flow net for a grid of complex potentials.
%
%   VARIABLES
%   X - x-values for grid (double vector)
%   Y - y-values for grid (double vector)
%   Grid - grid of complex potential (complex matrix)
%   lvs - number of contour levels (integer)
%   ln_wdth - line width (double)
%
%   LATEST UPDATE
%   2021-10-26
%
%   AUTHOR
%   Erik Toller,
%   Department of Earth Sciences, Uppsala University, Sweden

% Find the max and min value for the stream lines
im_max=max(max(imag(Grid)));
im_min=min(min(imag(Grid)));
% Calculate the difference between max an min
D=im_max-im_min;
% Calcualte the step size for the given number of levels
del=D/lvs;

% Find the max and min value for the stream lines
re_max=max(max(real(Grid)));
re_min=min(min(real(Grid)));
% Calculate the difference between max an min
D=re_max-re_min;
% Calcualte the number of levels using the same step size as for the
% stream lines
lvsr=round(D/del);

% Export the levels
lvs_im = linspace(im_min,im_max,lvs);
lvs_re = linspace(re_min,re_max,lvsr);

% Contour stream lines and equipotentials
M1 = contour(X, Y,real(Grid),lvs_re,'r','LineWidth',ln_wdth);
M2 = contour(X, Y,imag(Grid),lvs_im,'b','LineWidth',ln_wdth);
end