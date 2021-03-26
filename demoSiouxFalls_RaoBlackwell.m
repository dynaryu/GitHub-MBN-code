%{
Sioux Falls network connectivity analysis using inexhaustive CPM
whose rules are stochastically searched

P(S=0) & P(C*=0|S=0,I*=1)
%}
clear; close all;
import mbn.*
%% Generate samples
load SiouxFalls_BN % from 'SiouxFalls_BN.m'
load SiouxFalls_RDA sNode tNode
pickidx = 38;

nSample = 10e4; % # of samples
Px1 = .95; % P(c1) 
rng(777)
P.x = ( rand(nSample,nComp)<Px1 );
P.Q = Px1.^P.x .* (1-Px1).^(1-P.x);
P.Q = exp(sum(log(P.Q),2));

nNode = max(pairs(:));

%% Rao-Blackwell [1] P(S,i1)
D.w1 = zeros( nSample,1 ); % for P(S=0) 
D.w2 = zeros( nSample,1 ); % for P(C=0|S=0,I*=0)
D.mu1 = zeros( nSample,1 ); % for P(S=0) 
D.mu2 = zeros( nSample,1 ); % for P(C=0|S=0,I*=0)
Sfail = zeros(nSample,1); 

for ii = 1:nSample
    c_ = P.x(ii,:);
    
%     Connectivity check
    A = eye(nNode);
    for jj = 1:nComp
        if c_(jj) == 1
            A(pairs(jj,1),pairs(jj,2))=1;
        end
    end
    A = A^(nNode-1);
    
%     Construct cpm for current sample
    cpm_ = cpmSioux;  
    if A(sNode,tNode)
        cpm_{end}.C = [1 c_]; 
    else
        cpm_{end}.C = [0 c_]; 
        Sfail(ii) = 1;
    end
    cpm_{end}.p = 1;
    
%     Reduce-and-VE
    cpm1_ = conditioning( cpm_,var.C1-1+(1:nComp),c_ );
    cpm2_ = conditioning( cpm_,[var.C1-1+(1:nComp) var.S var.I1-1+pickidx],[c_ 0 1] );
    
    cpm1_( var.I1-1+(1:nComp) ) = [];
    cpm2_( var.I1-1+setdiff(1:nComp,pickidx) ) = [];
    
%     Inference approach [1]
%     cpm1_ = sumProductVE( cpm1_,var.D1-1+(1:nComp) );
%     cpm1_ = sum( cpm1_,[var.M var.L var.C1-1+(1:nComp)]);
%     cpm2_ = sumProductVE( cpm2_,var.D1-1+(1:nComp) );
%     cpm2_ = sum( cpm2_,[var.M var.L var.C1-1+setdiff(1:nComp,pickidx)]);

%     Inference approach [2]
    cpm1_ = sumProductVE( cpm1_,[var.D1-1+(1:nComp) var.M var.L var.C1-1+(1:nComp)] );    
    cpm2_ = sumProductVE( cpm2_,[var.D1-1+(1:nComp) var.M var.L var.C1-1+setdiff(1:nComp,pickidx)] );
    
%     Compute-Particle
    Pxe1 = cpm1_.p; D.w1(ii) = Pxe1/P.Q(ii); D.mu1(ii) = 1-cpm1_.C;
    
    if ~isempty( cpm2_.C )
        Pxe2 = cpm2_.p; D.w2(ii) = Pxe2/P.Q(ii);
        D.mu2(ii) = 1-cpm2_.C(1);
    end    
    
    if rem(ii,1e1) == 0        
        w1_ = D.w1(1:ii); w2_ = D.w2(1:ii); Sfail_ = Sfail(1:ii);
        mu1_ = D.mu1(1:ii); mu2_ = D.mu2(1:ii);
        
        disp( ['[**] ' num2str(ii) '-th Samples being done | # of failure event: ' num2str(sum(Sfail_))] );
        PS0 = mean( w1_.*mu1_ ); 
        PS0_var = mean(w1_.^2 .* mu1_.^2) - PS0^2;
        PS0_cov = sqrt( PS0_var/ii ) / PS0;
        disp([ 'Estimate of P(S0) = ' num2str( PS0 ) ' with c.o.v. ' num2str( PS0_cov )])
        
        f_ = P.x(1:ii,pickidx) == 0;
        PC0S0I1 = sum( w2_ .* mu2_ ) / sum(w2_);
        PC0S0I1_var = sum( w2_.^2 .* (mu2_-PC0S0I1).^2 ) / sum(w2_)^2;
        PC0S0I1_cov = sqrt(PC0S0I1_var/sum(w2_>0)) / PC0S0I1;
        disp([ 'Estimate of P(C0|S0,I1) = ' num2str( PC0S0I1 ) ' with c.o.v. ' num2str( PC0S0I1_cov )])     
        
        if rem(ii,10e3)==0
            eval( ['save demoSiouxFalls_RaoBlackwell_v5_' num2str(ii/10e3)] );
        end
        
    end
end

save demoSiouxFalls_RaoBlackwell

disp([ '95% Confidence interval of P(S0): [ ' num2str( PS0*(1-1.96*PS0_cov) ) ...
    ' ' num2str( PS0*(1+1.96*PS0_cov) ) ' ]'])
disp([ '95% Confidence interval of P(C0|S0,I1): [ ' num2str( PC0S0I1*(1-1.96*PC0S0I1_cov) ) ...
    ' ' num2str( PC0S0I1*(1+1.96*PC0S0I1_cov) ) ' ]'])
