function h = makeMyFigure(W,H)
% makeMyFigure(W,H) makes a figure that is W cm x H cm

% if one argument, use golden ratio.
if nargin == 1
    H = W;
    W = W*(1+sqrt(5))/2;
end
h = figure();
set(gcf, 'PaperPositionMode','manual', 'PaperUnits','centimeters', ...
            'PaperSize', [W H], 'PaperPosition',[0 0 W H], ...
            'Units','centimeters','Position',[0 0 W H]);