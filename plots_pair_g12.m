
%%
% clims = [0 0.5364];

%% 
data=readmatrix('pair_TG.txt');
xrange = data(:,1);
yrange = data(:,2);
fun = data(:,3:end);
clims = [0 max(max(fun))];
figure(101)
imagesc(xrange,yrange,fun,clims)
pbaspect([1 1 1])
set(gca,'FontSize',32)
set(gca, 'FontName', 'Times New Roman')
set(gca,'YDir','normal')
colorbar
% set(gca,'XTick',[], 'YTick', [])
xlim([-3.2 3.2])
ylim([-3.2 3.2])

%%
data=readmatrix('pair_g121_W26.txt');
xrange = data(:,1);
yrange = data(:,2);
fun = data(:,3:end);
clims = [0 max(max(fun))];
figure(102)
imagesc(xrange,yrange,fun,clims)
pbaspect([1 1 1])
set(gca,'FontSize',32)
set(gca, 'FontName', 'Times New Roman')
set(gca,'YDir','normal')
colorbar
set(gca,'XTick',[], 'YTick', [])
xlim([-1.6 1.6])
ylim([-1.6 1.6])

%%
data=readmatrix('pair_g1205_W26.txt');
xrange = data(:,1);
yrange = data(:,2);
fun = data(:,3:end);
clims = [0 max(max(fun))];
figure(103)
imagesc(xrange,yrange,fun,clims)
pbaspect([1 1 1])
set(gca,'FontSize',32)
set(gca, 'FontName', 'Times New Roman')
set(gca,'YDir','normal')
colorbar
set(gca,'XTick',[], 'YTick', [])
xlim([-1.6 1.6])
ylim([-1.6 1.6])

%%
data=readmatrix('pair_g122_W26.txt');
xrange = data(:,1);
yrange = data(:,2);
fun = data(:,3:end);
clims = [0 max(max(fun))];
figure(114)
imagesc(xrange,yrange,fun,clims)
pbaspect([1 1 1])
set(gca,'FontSize',32)
set(gca, 'FontName', 'Times New Roman')
set(gca,'YDir','normal')
colorbar
set(gca,'XTick',[], 'YTick', [])
xlim([-1.6 1.6])
ylim([-1.6 1.6])

%%
data=readmatrix('pair_g123_W26.txt');
xrange = data(:,1);
yrange = data(:,2);
fun = data(:,3:end);
clims = [0 max(max(fun))];
figure(104)
imagesc(xrange,yrange,fun,clims)
pbaspect([1 1 1])
set(gca,'FontSize',32)
set(gca, 'FontName', 'Times New Roman')
set(gca,'YDir','normal')
colorbar
set(gca,'XTick',[], 'YTick', [])
xlim([-1.6 1.6])
ylim([-1.6 1.6])

%%
data=readmatrix('pair_g121_W30.txt');
xrange = data(:,1);
yrange = data(:,2);
fun = data(:,3:end);
clims = [0 max(max(fun))];
figure(105)
imagesc(xrange,yrange,fun,clims)
pbaspect([1 1 1])
set(gca,'FontSize',28)
set(gca, 'FontName', 'Times New Roman')
set(gca,'YDir','normal')
colorbar
set(gca,'XTick',[], 'YTick', [])
xlim([-1.6 1.6])
ylim([-1.6 1.6])

%%
data=readmatrix('pair_g121_W20.txt');
xrange = data(:,1);
yrange = data(:,2);
fun = data(:,3:end);
clims = [0 max(max(fun))];
figure(106)
imagesc(xrange,yrange,fun,clims)
pbaspect([1 1 1])
set(gca,'FontSize',32)
set(gca, 'FontName', 'Times New Roman')
set(gca,'YDir','normal')
colorbar
set(gca,'XTick',[], 'YTick', [])
xlim([-1.6 1.6])
ylim([-1.6 1.6])

%%
data=readmatrix('pair_g121_W10.txt');
xrange = data(:,1);
yrange = data(:,2);
fun = data(:,3:end);
clims = [0 max(max(fun))];
figure(107)
imagesc(xrange,yrange,fun,clims)
pbaspect([1 1 1])
set(gca,'FontSize',32)
set(gca, 'FontName', 'Times New Roman')
set(gca,'YDir','normal')
colorbar
% set(gca,'XTick',[], 'YTick', [])
% xticks([-1.5 -1.0 -0.5 0 0.5 1.0 1.5])
xlim([-1.6 1.6])
ylim([-1.6 1.6])

%%









