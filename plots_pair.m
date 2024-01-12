
%%
% clims = [4 18];

%% 
data=readmatrix('pair_TG.txt');
xrange = data(:,1);
yrange = data(:,2);
fun = data(:,3:end);
figure(101)
imagesc(xrange,yrange,fun)
pbaspect([1 1 1])
set(gca,'FontSize',28)
set(gca, 'FontName', 'Times New Roman')
set(gca,'YDir','normal')
% colorbar
set(gca,'XTick',[], 'YTick', [])

% saveas(figure(101),'pair_TG.fig')
% saveas(figure(101),'pair_TG.png')

%%
data=readmatrix('pair_g121_W26.txt');
xrange = data(:,1);
yrange = data(:,2);
fun = data(:,3:end);
figure(102)
imagesc(xrange,yrange,fun)
pbaspect([1 1 1])
set(gca,'FontSize',28)
set(gca, 'FontName', 'Times New Roman')
set(gca,'YDir','normal')
% colorbar
set(gca,'XTick',[], 'YTick', [])

%%
data=readmatrix('pair_g1205_W26.txt');
xrange = data(:,1);
yrange = data(:,2);
fun = data(:,3:end);
figure(103)
imagesc(xrange,yrange,fun)
pbaspect([1 1 1])
set(gca,'FontSize',28)
set(gca, 'FontName', 'Times New Roman')
set(gca,'YDir','normal')
% colorbar
set(gca,'XTick',[], 'YTick', [])

%%
data=readmatrix('pair_g123_W26.txt');
xrange = data(:,1);
yrange = data(:,2);
fun = data(:,3:end);
figure(104)
imagesc(xrange,yrange,fun)
pbaspect([1 1 1])
set(gca,'FontSize',28)
set(gca, 'FontName', 'Times New Roman')
set(gca,'YDir','normal')
% colorbar
set(gca,'XTick',[], 'YTick', [])


%%









