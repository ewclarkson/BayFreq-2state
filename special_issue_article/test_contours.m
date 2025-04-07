%% test iso-contour plots for gaussian random data

Ndata = 1E5; % number of data points in each direction
x = randn(Ndata,1);
y = randn(Ndata,1);
[bc,Xedges,Yedges] = histcounts2(x,y,'Normalization','probability');
bar3(bc,1)

Pmax = max(bc,[],'all'); % define maximum
lv1 = Pmax*exp(-1/2);
lv2 = Pmax*exp(-2^2/2); 
lv3 = Pmax*exp(-3^2/2);
% lv1 = 4*10^(-3);
% lv2 = 1*10^(-3);
% lv3 = 1*10^(-4);

contourf(bc, [lv1,lv1],'EdgeColor','none','FaceColor',[0 0.3 1],'FaceAlpha',0.3)
hold on
contourf(bc, [lv2,lv2],'EdgeColor','none','FaceColor',[0 0.3 1],'FaceAlpha',0.3)
hold on
contourf(bc, [lv3,lv3],'EdgeColor','none','FaceColor',[0 0.3 1],'FaceAlpha',0.3)













