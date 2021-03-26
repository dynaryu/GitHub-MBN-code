function [pathNode,rel,pathLink] = mDijkstra( Adj,Anumber,Prob,sNode,tNode )
%MDIJKSTRA Get the most reliable path for sRDA (Lim and Song 2012)
% [path,rel] = mDijkstra( Adj,Prob,sNode,tNode )
%{
Links are subject to failure (NOT nodes.)
%}

%{
Adj: adjacent matrix
Anumbered: N x N matrix w/ links' numbering (N: # of nodes)
Prob: N x N matrix where only the non-zero elements in Adj has non-zero
probability of being in operation.
sNode: source node
tNode: terminal node
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
%}

%%
nNode = size(Adj,1);
S = sNode; % parmanantly labeled
S_ = setdiff(1:nNode,S); % temporarily labeled
r = zeros(1,nNode); r(sNode) = 1;
pred = zeros(1,nNode);
for ii = find( Adj(sNode,:) )
    r(ii) = Prob(sNode,ii);
    pred(ii) = sNode;
end
while length(S) < nNode
    [~,ii] = sort( r,'descend' );
    ii = setdiff( ii,S,'stable' ); ii = ii(1);

    S = union(S,ii); S_ = setdiff(S_,ii);
    
    jj_ = find( Adj(ii,:) );
    for jj = jj_
        if r(jj) < r(ii)*Prob(ii,jj)
            r(jj) = r(ii)*Prob(ii,jj);
            pred(jj) = ii;
        end
    end
end
if ~r(tNode)
    pathNode = []; rel = 0; pathLink = [];
else
    jj = tNode; pathNode = jj; rel = 1; pathLink = [];
    while jj ~= sNode
        pathNode = [pred(jj) pathNode];
        rel = rel*Prob(pred(jj),jj);
        pathLink = [Anumber(pathNode(1),pathNode(2)) pathLink];
        jj = pred(jj);
    end
end