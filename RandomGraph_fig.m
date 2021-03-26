clear; close all;

load demoRandomGraph
%% plot1: # of elements required for CPT
Fsz_ppt = 16; Fsz_ppr = 13;
Fsz_ppt_label = 14; Fsz_ppr_label = 11;
LW = 1.5;
nLINK_ = nLINK(1:length(nGRAPH));

figure;
plot( nLINK_,nGRAPH,'*--','Linewidth',LW )
grid on
ax = gca;
ax.XAxis.FontSize = Fsz_ppt_label;
ax.YAxis.FontSize = Fsz_ppt_label;
xlabel( 'Number of links','Fontsize',Fsz_ppt,'FontName','times new roman' )
ylabel( 'Number of identified subsets','Fontsize',Fsz_ppt,'FontName','times new roman' )
saveas(gcf,'figure/CutLinkSets_GraphNo.emf')
saveas(gcf,'figure/CutLinkSets_GraphNo.pdf')

figure;
semilogy( nLINK_,nElem_MSR,'s--','Linewidth',LW )
hold on
semilogy(nLINK_,2.^nLINK_.*(nLINK_+1),'*--','Linewidth',LW )
grid on
ax = gca;
ax.XAxis.FontSize = Fsz_ppt_label;
ax.YAxis.FontSize = Fsz_ppt_label;
xlabel( 'Number of links','Fontsize',Fsz_ppt,'FontName','times new roman' )
ylabel( 'Number of parameters','Fontsize',Fsz_ppt,'FontName','times new roman' )
% axis([nLINK(1)-2 nLINK_(end)+2 1 10^20])
legend({'MBN','CPT of naive formulation'},'Fontsize',Fsz_ppt,'Location','NorthWest','FontName','times new roman')
saveas(gcf,'figure/CutLinkSets_CPMsize.emf')
saveas(gcf,'figure/CutLinkSets_CPMsize.pdf')


%% plot2: Graph topology

s = []; t = [];
Anumber_ = Anumber;
for ii = 1:nLINK_(end)
    [s_,t_] = find(Anumber_==ii);
    s = [s; s_];
    t = [t; t_];
end

figure('rend','painters','pos',[10 10 700 600])
Graph = digraph(s,t);
nnode_ = max( [max(s),max(t)] );
for ii = 1:nnode_
    nlabel{ii,1} = num2str(ii);
end
nlabel{sNode} = 's(1)'; nlabel{tNode} = 't(11)';
h = plot(Graph,'Layout','circle','NodeLabel',{},'LineWidth',3,'ArrowSize',15);
markloc_ = 1.1*[cos( 2*pi/nnode_*(1:nnode_)' ) sin( 2*pi/nnode_*(1:nnode_)' )];
markloc_(8:14,:) = markloc_(8:14,:)*1.15/1.1;
markloc_(10,:) = markloc_(10,:)*1.2/1.15; 
markloc_(11,:) = markloc_(11,:)*1.305/1.15; 
for ii=1:nnode_
   text(markloc_(ii,1),markloc_(ii,2),nlabel{ii},'fontsize',Fsz_ppt,'FontWeight','bold','FontName','times new roman');
end
ax = gca;
ax.Visible = 'off';
saveas(gcf,'figure/CutLinkSets_network.emf')
saveas(gcf,'figure/CutLinkSets_network.pdf')

%% Plot 3: Inference
figure;
plot( nLINK_,PS0,'*--','Linewidth',LW )
hold on
plot( nLINK_,PC0S0,'*--','Linewidth',LW )
grid on
ax = gca;
ax.XAxis.FontSize = Fsz_ppt_label;
ax.YAxis.FontSize = Fsz_ppt_label;
xlabel( 'Number of links','Fontsize',Fsz_ppt,'FontName','times new roman' )
ylabel( 'Probability','Fontsize',Fsz_ppt,'FontName','times new roman' )
legend({'P(X_{N+1}= 0)','P(X_{22}= 0 | X_{N+1}= 0)'},'Fontsize',Fsz_ppt-2,...
    'Location','best','FontName','times new roman')
saveas(gcf,'figure/InfResult.emf')
saveas(gcf,'figure/InfResult.pdf')