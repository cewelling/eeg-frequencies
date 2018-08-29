% low dominant segments 
lSegs_low = lowSegList(strcmp(labelVect, 'lo'),:); 
lSegs_high = highSegList(strcmp(labelVect, 'lo'),:);

% high dominant segments
hSegs_low = lowSegList(strcmp(labelVect, 'hi'),:);
hSegs_high = highSegList(strcmp(labelVect, 'hi'),:);

c=[ones(1,size(lSegs_low, 1)) 2*ones(1,size(lSegs_high,1))];
x=1:size(lowSegList, 2);
y=[lSegs_low;lSegs_high];

clear g
g(1,1)=gramm('x',x,'y',y,'color',c);
g(1,2)=gramm('x',x,'y',y,'color',c);

c=[ones(1,size(hSegs_low, 1)) 2*ones(1,size(hSegs_high,1))];
x=1:size(highSegList, 2);
y=[hSegs_low;hSegs_high];

g(2,1) = gramm('x',x,'y',y,'color',c);
g(2,2) = gramm('x',x,'y',y,'color',c);

g(1,1).geom_line();
g(1,1).set_title('Individual Segments');
%g(1,1).set_color_options('map','matlab'); %It is also possible to provide a custom
% colormap by providing a N-by-3 matrix (columns are R,G,B).

g(1,2).stat_summary();
g(1,2).set_title('Group Averaged Segments');
%g(1,2).set_color_options('map','matlab');

g(2,1).geom_line();
g(2,1).set_title('Individual Segments');
%g(1,1).set_color_options('map','matlab'); %It is also possible to provide a custom
% colormap by providing a N-by-3 matrix (columns are R,G,B).

g(2,2).stat_summary();
g(2,2).set_title('Group Averaged Segments');
%g(1,2).set_color_options('map','matlab');

g.set_title('Rivalry Segments for classification analysis');

figure('Position',[100 100 800 550]);
g.draw();
