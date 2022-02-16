%% Full Network
DG = distances(GD);
[diam, eccs] = diams(DG);
[eff, hets] = effs(DG);
num = numnodes(GD);

GD2 = GD.rmedge(1,6);
GD2 = GD2.rmedge(6,1);
DG2 = distances(GD2);
[diam2, ~] = diams(DG2);
UGC2 = UGD.rmedge(1,6);
num2 = numnodes(GD2);

GD3 = GD2.rmedge(1,4);
GD3 = GD3.rmedge(4,1);
GD3 = GD3.rmedge(3,1);
GD3 = GD3.rmedge(1,3);
DG3 = distances(GD3);
[diam3, ~] = diams(DG3);
UGC3 = UGC2.rmedge(1,4);
UGC3 = UGC3.rmedge(1,3);
num3 = numnodes(GD3);
% cmap = flipud(autumn);
weights = ones(length(UGD.Edges.Weight),1);
weights = weights * 4;

%% Plotting    
%set(gcf, 'Position',  posize);
%set(gcf,'color',color);

lwidth = 3*UGD.Edges.Weight/max(UGD.Edges.Weight);
lwidth2 = 3*UGC2.Edges.Weight/max(UGC2.Edges.Weight);
lwidth3 = 3*UGC3.Edges.Weight/max(UGC3.Edges.Weight);

figure(1)

subplot(1,3,1) %%%%%%%%%%%%%%
p=plot(UGC,'XData',xlocation,'YData',ylocation,...
'EdgeLabel',UGD.Edges.Weight,...
'NodeLabel', 1:num,...
'LineWidth',lwidth);...
title('Nodal Diameters');

UGD.Nodes.NodeColors = diam;
p.NodeCData = UGD.Nodes.NodeColors;

nl = p.NodeLabel;
p.NodeLabel = '';
xd = get(p, 'XData');
yd = get(p, 'YData');
text(xd, yd, nl, 'FontSize',10,...
'HorizontalAlignment','right', 'VerticalAlignment','bottom')
colormap autumn; colorbar
    
subplot(1,3,2) %%%%%%%%%%%%%%

p=plot(UGC2,'XData',xlocation,'YData',ylocation,...
'EdgeLabel',UGC2.Edges.Weight,...
'NodeLabel', 1:num2,...
'LineWidth',lwidth2);...
title('Nodal Diameters');

UGD.Nodes.NodeColors = diam2;
p.NodeCData = UGD.Nodes.NodeColors;

nl = p.NodeLabel;
p.NodeLabel = '';
xd = get(p, 'XData');
yd = get(p, 'YData');
text(xd, yd, nl, 'FontSize',10,...
'HorizontalAlignment','right', 'VerticalAlignment','bottom')
colormap autumn; colorbar


subplot(1,3,3) %%%%%%%%%%%%%%

p=plot(UGC3,'XData',xlocation,'YData',ylocation,...
'EdgeLabel',UGC3.Edges.Weight,...
'NodeLabel', 1:num3,...
'LineWidth',lwidth3);...
title('Nodal Diameters');

UGD.Nodes.NodeColors = diam3;
p.NodeCData = UGD.Nodes.NodeColors;

nl = p.NodeLabel;
p.NodeLabel = '';
xd = get(p, 'XData');
yd = get(p, 'YData');
text(xd, yd, nl, 'FontSize',10,...
'HorizontalAlignment','right', 'VerticalAlignment','bottom')
colormap autumn; colorbar



