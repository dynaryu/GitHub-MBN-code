%{
Examines the # of parameters required for the MBN to describe the
connectivity of a given graph
%}
import mbn.*
clear; close all;
%% Generate random-graph
nNODE = 10+(1:14);
nLINK = floor(nNODE.*(nNODE-1)*.5 * .2*.2)*5;
nLINK(1:8) = 22:29;
nodeLinkRat = nLINK./(nNODE.\(nNODE-1)*.5);

rng(1)
nNode_ = nNODE(1); nLink_ = nLINK(1);
sNode = 1; tNode = nNode_;
Adj = zeros( nNode_ ); Anumber = zeros( nNode_ ); Prob = zeros( nNode_ );
idx_ = sub2ind( size(Adj),1:nNode_,1:nNode_ );
idx_ = setdiff( 1:nNode_^2,idx_ );
idxpick_ = randsample(idx_,nLink_);
idxpick_ = sort(idxpick_);
Adj(idxpick_) = 1; Anumber(idxpick_) = 1:nLink_; 
Prob(Adj>0) = rand( nLink_,1 );
ProbLink = Prob(Prob>0);

[G,sysR,sysF,err] = sRDA(Adj,Anumber,Prob,ProbLink,sNode,tNode,0,1,1e5);
nGRAPH = length(G.nEvent);
nElem_MBN = length(G.nEvent)*(nLink_+2);
SYSR = sysR(end);
nEVENT = zeros(length(nLINK),1);
nEVENT(1) = sum(G.nEvent);
nEVENTexact = 2.^(nLINK+1)';

cpick_ = 22;
PS0 = zeros( length(nLINK),1 ); % P(S0)
PC0S0 = zeros( length(nLINK),1 ); % P(Cpick=0|S0)
% BN Inference
for jj = 1:nLink_
    cpmC{jj,1} = cpmjoint( jj,[0;1],[0.1;0.9],{'Fail','Survive'} );
end
SC_ = [];
for ii = 1:length(G.S)
    SCtmp_ = -ones(1,nLink_+1);
    if ~isempty(G.PATHLINK{ii})
        SCtmp_(1) = 1;
    else
        SCtmp_(1) = 0;
    end
    SCtmp_( G.PATHLINK{ii}+1 ) = 1;
    SCtmp_( G.S{ii} ) = 1;
    SCtmp_( G.S_{ii} ) = 0;
    
    SC_ = [SC_; SCtmp_];
end
cpmS{1} = cpmcond( nLink_ + 1,1:nLink_,SC_,ones( length(G.S),1 ) );
cpmSmarginal = sumProductVE( [cpmC; cpmS],1:nLink_ );
cpmSC1 = sumProductVE( conditioning([cpmC; cpmS],nLink_+1,0),setdiff(1:nLink_,cpick_) );
PS0(1) = cpmSmarginal.p(1);
PC0S0(1) = sum(cpmSC1.p(1:2)) / sum( cpmSC1.p );


for kk = 2:length(nNODE)
    nNodep_ = nNode_; nLinkp_ = nLink_;
    nNode_ = nNODE(kk); nLink_ = nLINK(kk);
    nNodeadd_ = nNode_-nNodep_; nLinkadd_ = nLink_-nLinkp_;
    Adj = [Adj zeros(nNodep_,1)]; Adj = [Adj; zeros(1,nNode_)];
    Anumber = [Anumber zeros(nNodep_,1)]; Anumber = [Anumber; zeros(1,nNode_)];
    Prob = [Prob zeros(nNodep_,1)]; Prob = [Prob; zeros(1,nNode_)];
    
    idx_ = setdiff( find(~Adj),sub2ind(size(Adj),1:nNode_,1:nNode_ ) );

    idxpick_ = randsample(idx_,nLinkadd_);
    Adj(idxpick_) = 1; Anumber(idxpick_) = nLinkp_+1:nLink_; 
    ProbLink = [ProbLink; rand( nLinkadd_,1 )];
    Prob(idxpick_) = ProbLink(nLinkp_+1:end);
    
    [G,sysR,sysF,err] = sRDA(Adj,Anumber,Prob,ProbLink,sNode,tNode,0,1,2e4);
    if length(err) == 2e4+1
        break;
    end
    nElem_MSR(kk,1) = length(G.nEvent)*(nLink_+2); % elements of c vector + p vector
    SYSR(kk,1) = sysR(end);
    nGRAPH(kk,1) = length(G.nEvent);
    nEVENT(kk,1) = sum( G.nEvent);
    
    % BN Inference
    for jj = (nLinkp_+1):nLink_
        cpmC{jj,1} = cpmjoint( jj,[0;1],[0.1;0.9],{'Fail','Survive'} );
    end
    SC_ = [];
    for ii = 1:length(G.S)
        SCtmp_ = -ones(1,nLink_+1);
        if ~isempty(G.PATHLINK{ii})
            SCtmp_(1) = 1;
        else
            SCtmp_(1) = 0;
        end
        SCtmp_( G.PATHLINK{ii}+1 ) = 1;
        SCtmp_( G.S{ii}+1 ) = 1;
        SCtmp_( G.S_{ii}+1 ) = 0;

        SC_ = [SC_; SCtmp_];
    end
    cpmS{1} = cpmcond( nLink_ + 1,1:nLink_,SC_,ones( length(G.S),1 ) );
    cpmSmarginal = sumProductVE( [cpmC; cpmS],1:nLink_ );
    cpmSC1 = sumProductVE( conditioning([cpmC; cpmS],nLink_+1,0),setdiff(1:nLink_,cpick_) );
    PS0(kk) = cpmSmarginal.p(1);
    PC0S0(kk) = sum( cpmSC1.p(1:2) ) / sum( cpmSC1.p );   

end
   
% save demoRandomGraph

%% Figures
% run( 'RandomGraph_fig.m' )
