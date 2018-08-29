%little figure

box off;
set(gca,'TickDir','out');
colormap(parula(5))


xh = get(gca,'xlabel') % handle to the label object
p = get(xh,'position') % get the current position property
p(2) = 1.5*p(2) ;        % double the distance, 
set(xh,'position',p)   % set the new position

yh = get(gca,'ylabel') % handle to the label object
p = get(yh,'position') % get the current position property
p(1) = p(1) - 0.35*p(1) ;        % double the distance, 
set(yh,'position',p)   % set the new position
