% function demoSiouxFalls_RDA
%{
Sioux Falls network connectivity analysis using inexhaustive CPM
whose rules are deterministically searched by RDA

P(C*=0|S=0,I*=1)
%}
clear; close all;
import mbn.*
%% load data 
load SiouxFalls_BN % from 'SiouxFalls_BN.m'
pickidx = 38;

nSample = length(cpmSioux{var.S}.p); % # of samples
SCC = cpmSioux{var.S}.C(:,2:end);
SC1 = cpmSioux{var.S}.C(:,1);

%% [1] Eliminate I's - {I*}
cpmSioux2 = cpmSioux;
cpmSioux2{var.D1-1+pickidx} = product( cpmSioux{var.D1-1+pickidx},cpmSioux{var.I1-1+pickidx} );
cpmSioux2(var.I1-1+(1:nComp)) = [];

%% [2] Eliminate D's
cpmSioux3 = cpmSioux2;
for ii = 1:nComp
    cpmSioux3{var.I1-1+ii} = product( cpmSioux2{var.D1-1+ii},cpmSioux2{var.I1-1+ii} );
    cpmSioux3{var.I1-1+ii} = sum( cpmSioux3{var.I1-1+ii},var.D1-1+ii);
end
cpmSioux3(var.D1-1+(1:nComp)) = [];

%% [3] Get P(S,C*,I*)
scope5 = [var.S var.C1-1+pickidx var.I1-1+pickidx];

scope4 = [var.S var.C1-1+pickidx var.I1-1+pickidx var.M var.L]; % For conditioning on M,L
cpmML = product( cpmSioux{var.M},cpmSioux{var.L} );

C5 = []; p5 = [];
for ii = 1:size(cpmML.C,1)
       
    ml_ = cpmML.C(ii,:);
    cpm_ = conditioning( cpmSioux3,[var.M var.L],ml_);
    
    p4_ = zeros( nSample,1) + log( cpmML.p(ii) );% P(M)P(L)
    for jj = setdiff(1:nComp,pickidx)
        p4_ = p4_ + (SCC(:,jj)==0) * log( cpm_{2+jj}.p(1) ) + (SCC(:,jj)==1) * log( cpm_{2+jj}.p(2) );
    end
    p4_ = exp(p4_);    
    
    cpm2_ = cpmcond( var.S,var.C1-1+pickidx,[SC1 SCC(:,pickidx)],p4_ );
    cpm2_ = product( cpm2_,cpm_{pickidx+2} );
    p4 = []; C4 = []; C2_ = cpm2_.C; p2_ = cpm2_.p;
    while ~isempty( p2_ )
        c2_ = C2_(1,:);
        idx_ = ismember( C2_,c2_,'rows' );
        
        p4 = [p4; sum( p2_(idx_) )];
        C4 = [C4; c2_];
        
        p2_(idx_) = [];
        C2_(idx_,:) = [];
    end
    
    cpmSioux4{ii} = cpmjoint( [cpm2_.scope cpm2_.scopep],C4,p4 );
    
    C5 = [C5; C4]; p5 = [p5; p4];
end

cpmSioux5 = cpmjoint( cpmSioux4{1}.scope,C5,p5 );
cpmSioux5 = sum( cpmSioux5,[var.M var.L] );

cpmSioux6 = sum( cpmSioux5,var.C1-1+pickidx,1 );
    
save demoSiouxFalls_RDA

%% Probability
idx1_ = cpmSioux5.C(:,2) == 0; % S=0
PS0_low = sum( cpmSioux5.p(idx1_ ) );
PS0_up = 1-sum( cpmSioux5.p(~idx1_ ) );
disp([ 'P(S0) = [ ' num2str( PS0_low ) ' ' num2str( PS0_up ) ' ]' ])

idx2_ = ismember( cpmSioux5.C(:,2:3),[0 1],'rows' ); % (s0,i1)
PS0I1_low = sum( cpmSioux5.p(idx2_ ) );
PS0I1_up = 1-sum( cpmSioux5.p(~idx2_ ) );
disp([ 'P(S0,I1) = [ ' num2str( PS0I1_low ) ' ' num2str( PS0I1_up ) ' ]' ])

idx3_ = ismember( cpmSioux5.C,[0 0 1],'rows' ); % (c0,s0,i1)
PC1S0I1_low = sum( cpmSioux5.p(idx3_ ) );
PC1S0I1_up = 1-sum( cpmSioux5.p(~idx3_ ) );
disp( [ 'P(C0|S0,I1) = [ ' num2str( PC1S0I1_low/PS0I1_up ) ' ' num2str( PC1S0I1_up/PS0I1_low ) ' ]' ] )

disp( ['P(C0) = ' num2str( cpmSioux6.p(2) ) ] )
