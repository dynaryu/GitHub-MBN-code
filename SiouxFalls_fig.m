% figure for Sioux Falls network example

Fsz_ppt = 16; Fsz_ppr = 13;
Fsz_ppt_label = 14; Fsz_ppr_label = 11;
LW = 1.5;

figure('rend','painters','pos',[10 10 550 500])
h1 = plot(EQloc(:,1), EQloc(:,2),'ro','MarkerFaceColor','r','Markersize',6);
hold on
Graph = digraph(SFnet.arcinfo(:,1),SFnet.arcinfo(:,2));
h=plot(Graph,'XData',SFnet.locnode(:,1),'YData',SFnet.locnode(:,2),'LineWidth',LW,'EdgeLabel',1:length(SFnet.arcinfo));
h.NodeLabel = '';
el = h.EdgeLabel;
h.EdgeLabel = '';
xd = get(h, 'XData'); xd__ = xd([13 2]); xd__(1) = xd__(1)-1; xd__(2) = xd__(2)-.5; 
yd = get(h, 'YData'); yd__ = yd([13 2]); yd__(1) = yd__(1)-.3; yd__(2) = yd__(2)+.5; 
text(xd__, yd__, {'Source','Terminal'}, 'FontSize',Fsz_ppr_label,'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','middle','FontName','times new roman')

xd_ = xd( SFnet.arcinfo(:,1) ) + xd( SFnet.arcinfo(:,2) ); xd_ = .5*xd_;
yd_ = yd( SFnet.arcinfo(:,1) ) + yd( SFnet.arcinfo(:,2) ); yd_ = .5*yd_;
uplinks = [3 8 11 15 21 17 33 27 48 55 44 57 70 63 56 74 62 66];
downlinks = [1 6 9 12 24 20 36 32 29 50 30 41 45 60 72 68 39 75 64];
leftlinks = [2 7 37 10 34 42 73 13 25 28 46 69 16 22 49 53 59 18 56 4];
rightlinks = [5 35 38 31 40 71 76 23 26 43 67 65 14 19 47 52 58 61 60 54];
leftlinks2 = [37 10 34 42 73 41 70 46 72 66 69 68 51 49 12 16 21 24 48 22 49 53 39 56 18 55 20 18];
leftlinks3 = [33 36 74 11 27 32 44 75 64 68 45 57 51 15 17 50 28 59 13 25];
yd_(uplinks) = yd_(uplinks)+.5;
yd_(downlinks) = yd_(downlinks)-.5;
xd_(leftlinks) = xd_(leftlinks)-.5;
xd_(rightlinks) = xd_(rightlinks)+.2;
xd_(leftlinks2) = xd_(leftlinks2)-.2;
xd_(leftlinks3) = xd_(leftlinks3)-.2;
yd_(29) = yd_(29)+.2;
yd_(68) = yd_(68)+.2;
xd_(68) = xd_(68)-.2;
text(xd_, yd_, el, 'FontSize',Fsz_ppr_label,'FontWeight','bold', ...
    'HorizontalAlignment','left', 'VerticalAlignment','middle','FontName','times new roman')

xt = get(gca, 'XTick');
set(gca, 'FontSize', 12,'FontName','times new roman')
xlabel( 'x-direction (km)','Fontsize',Fsz_ppr,'FontName','times new roman' )
ylabel( 'y-direction (km)','Fontsize',Fsz_ppr,'FontName','times new roman' )
legend(h1,{'Epicenter locations'},'Fontsize',Fsz_ppr,'FontName','times new roman','Location','SouthEast')
saveas(gcf,'figure/SiouxFallsNetwork.pdf')
saveas(gcf,'figure/SiouxFallsNetwork.emf')