function [G,sysR,sysF,err] = sRDA(Adj,Anumber,Prob,ProbLink,sNode,tNode,tol,search,graphnum)
%SRDA: implementation of RDA (Li & He 2002) and sRDA (Lim and Song 2012)
% [G,sysR,sysF,err] = sRDA(Adj,Anumber,Prob,ProbLink,sNode,tNode,tol,<search>,<graphnum>)
%{
Links are subject to failure (NOT nodes.)
%}

%{
Adj: adjacent matrix
Anumbered: N x N matrix w/ links' numbering (N: # of nodes)
Prob: N x N matrix where only the non-zero elements in Adj has non-zero
probability of being in operation.
ProbLink: an array of probability of links according to link numbering
sNode: source node
tNode: terminal node
tol: allowed width for error bound (0: exact identification)
search: search mode -- 'breadth' or 1:breadth-first (preferable) / 'depth' or 2:depth-first
graphnum: the limit in # of graphs to be identified in case 'tol' cannot be acheived

G: structure w/ information on link sets and cut sets
    S: for each graph, the comp's in operation / S_: comp's in failure
    Adj: the adjacent matrix for each graph
    coef: coefficients of each graph
    PATHNODE: path of nodes for each link set / PATHLINK: path of links
sysR: an array of the lower bound of system reliability at each iteration
sysF: an array of the lower bound of system failure prob at each iteration
err : an array of the error width at each iteration 
%}

%{
e.g.,
pLink = 0.2; % probability of a link to be made
nNode = 10;
rng(1)
Adj = zeros(nNode);
Adj( rand(nNode)<pLink ) = 1; % Adjacent matrix
nLink = sum( Adj(:) );
Prob = Adj; % Probability of each link
Prob(Prob>0) = rand( nLink,1 );
sNode = 1; tNode = nNode;
[G,sysR,sysF,err] = sRDA(Adj,Anumber,Prob,ProbLink,sNode,tNode,0.01,1,1e4)
%}
import mbn.*
%%
if nargin<7
    search = 1;
elseif ischar(search)
    switch search
        case 'breadth'
            search = 1;
        case 'depth'
            search = 2;
    end
end
if nargin < 8
    graphnum = inf;
end        

%%
sysR = 0; sysF = 0;
Gidx = 1; % Current Graph being considered.

nLink = length(ProbLink);
G.S{Gidx} = []; G.S_{Gidx} = []; % S: comps in operation, S_: comps in failure
G.nEvent = []; % # of events considered.
G.Adj{Gidx} = Adj; G.coef(Gidx) = 1;
err = 1; nev_ = 0;

if search == 1

%     while ( err(end) > tol ) && ( length(err) < graphnum ) 
    while max( [( (tol>0) && ( err(end) > tol ) && (length(err) < graphnum+1) ), ...
            ( tol==0 ) && ( sum(G.nEvent) < 2^(nLink+1) ) && (length(err) < graphnum+1) ] )
        [pathNode rel pathLink] = mDijkstra( G.Adj{Gidx},Anumber,Prob,sNode,tNode );
        G.PATHNODE{Gidx,1} = pathNode;
        G.PATHLINK{Gidx,1} = pathLink;

        if ~isempty(pathLink)
            if isempty( intersect(pathLink,G.S{Gidx}) )
                sysr_ = log(rel) + log(G.coef(Gidx));
                nev_ = length(G.S{Gidx}) + length(G.S_{Gidx}) +length(pathLink);
            else                
                sysr_ = log(rel) - log( sum(ProbLink(intersect(pathLink,G.S{Gidx}))) ) + ...
                    log(G.coef(Gidx));
                nev_ = length(G.S{Gidx}) + length(G.S_{Gidx}) +length(pathLink) - ...
                    length( intersect(pathLink,G.S{Gidx}) );
            end
            sysR = [sysR sysR(end)+exp(sysr_)];
%             sysR = [sysR sysR(end)+rel*G.coef(Gidx)/prod(ProbLink(intersect(pathLink,G.S{Gidx})))];
            sysF = [sysF sysF(end)];
        else
            sysF = [sysF sysF(end)+G.coef(Gidx)];
            sysR = [sysR sysR(end)];
            nev_ = length(G.S{Gidx}) + length(G.S_{Gidx});
        end
        nev_ = nLink - nev_;
        nev_ = 2^(nev_+1); % for system event being both 0 and 1
        G.nEvent = [G.nEvent; nev_];

        pathLink = setdiff( pathLink,G.S{Gidx},'stable' );
        for ii = 1:length(pathLink)
            if ii == 1
                G.S{end+1,1} = G.S{Gidx};
            else
                G.S{end+1,1} = union(G.S{Gidx},pathLink(1:ii-1));
            end    
            G.S_{end+1,1} = union(G.S_{Gidx},pathLink(ii));
            adj_ = G.Adj{Gidx};
%             adj_( pathNode(ii),pathNode(ii+1) ) = 0;
            adj_ ( Anumber == pathLink(ii) ) = 0;
            
            G.Adj{end+1,1} = adj_;   
%             G.coef = [G.coef; prod(ProbLink(G.S{end}))*prod(1-ProbLink(G.S_{end}))];
            coef_ = sum( log(ProbLink(G.S{end})) ) + sum( log(1-ProbLink(G.S_{end})) );
            G.coef = [G.coef; exp(coef_)];
        end
        Gidx = Gidx+1;
        err = [err 1-sysF(end)-sysR(end)];
        
%         if sum(G.nEvent) == 2^(nLink+1); break; end
    end
    
else
    
    while ( err(end) > tol ) && (length(err) < graphnum)
        [pathNode rel pathLink] = mDijkstra( G.Adj{Gidx},Anumber,Prob,sNode,tNode );
        G.PATHNODE{Gidx,1} = pathNode;
        G.PATHLINK{Gidx,1} = pathLink;

        if ~isempty(pathLink)
            sysR = [sysR; sysR(end)+rel*G.coef(Gidx)/prod(ProbLink(intersect(pathLink,G.S{Gidx})))];
            sysF = [sysF; sysF(end)];

            pathLink = setdiff( pathLink,G.S{Gidx} );    
            if ~isempty(pathLink)        
                for ii = 1:length(pathLink)
                    if ii == 1
                        G.S = {G.S{1:Gidx} G.S{Gidx} G.S{Gidx+1:end}};
                    else
                        G.S = {G.S{1:Gidx+ii-1} union(G.S{Gidx},pathLink(1:ii-1)) G.S{Gidx+ii:end}};
                    end    
                    G.S_ = {G.S_{1:Gidx+ii-1} union(G.S_{Gidx},pathLink(ii)) G.S_{Gidx+ii:end}};
                    adj_ = G.Adj{Gidx};
                    adj_( pathNode(ii),pathNode(ii+1) ) = 0;
                    G.Adj = {G.Adj{1:Gidx+ii-1} adj_ G.Adj{Gidx+ii:end}};
%                     G.coef = [G.coef(1:Gidx+ii-1); prod(ProbLink(G.S{Gidx+ii}))*prod(1-ProbLink(G.S_{Gidx+ii}));...
%                         G.coef(Gidx+ii:end)];
                    coef_ = sum( log(ProbLink(G.S{end})) ) + sum( log(1-ProbLink(G.S_{end})) );
                    G.coef = [G.coef(1:Gidx+ii-1); exp(coef_);...
                        G.coef(Gidx+ii:end)];                    
                end
                Gidx = Gidx+1;   
            end
        else
                sysF = [sysF; sysF(end)+G.coef(Gidx)];
                sysR = [sysR; sysR(end)];
                Gidx = Gidx+1; 
        end

        err = [err; 1-sysF(end)-sysR(end)];
    end
end