%{
Sioux Falls network connectivity analysis by MCS

P(C*=0|S=0,I*=1)
%}
clear; close all;
import mbn.*
%% load data 
load SiouxFalls_BN % from 'SiouxFalls_BN.m'
load SiouxFalls_RDA sNode tNode
nNode= max(pairs(:));

pickidx = 38;

nVar = length(cpmSioux);
PS0 = []; PS0I1 = []; PC0S0I1=[];
nSample = 1e6;
covS0 = 1; covC0S0I1 = 1;

for ss = 1:nSample
    sampleIdx_ = zeros(nVar,1);
    cpmSioux_ = cpmSioux;

    % M, L
    sampleIdx_(1) = datasample( 1:5,1,'Weights',cpmSioux_{1}.p );
    sampleIdx_(2) = datasample( 1:10,1,'Weights',cpmSioux_{2}.p );
    cpmSioux_(var.C1-1+(1:nComp)) = conditioning( cpmSioux_(var.C1-1+(1:nComp)),[var.M var.L],sampleIdx_(1:2) );


    % D
    for ii = 1:nComp
        sampleIdx_(var.D1-1+ii) = (rand<cpmSioux_{var.D1-1+ii}.p(2));
        cpmSioux_{var.I1-1+ii} = conditioning( cpmSioux_{var.I1-1+ii},var.D1-1+ii,sampleIdx_(var.D1-1+ii) );
        cpmSioux_{var.C1-1+ii} = conditioning( cpmSioux_{var.C1-1+ii},var.D1-1+ii,sampleIdx_(var.D1-1+ii) );
    end

    % I, C
    for ii = 1:nComp
        sampleIdx_(var.I1-1+ii) = (rand<cpmSioux_{var.I1-1+ii}.p(2));
        sampleIdx_(var.C1-1+ii) = (rand<cpmSioux_{var.C1-1+ii}.p(2));
    end

    %     Connectivity check
    A = eye(nNode);
    for ii = 1:nComp
        if sampleIdx_(var.C1-1+ii) == 1
            A(pairs(ii,1),pairs(ii,2))=1;
        end
    end
    A = A^(nNode-1);
    if A(sNode,tNode)
        sampleIdx_(var.S) = 1;
    end

    PS0 = [PS0; 1-sampleIdx_(var.S)];
    PS0I1 = [PS0I1; PS0(end)*sampleIdx_(var.I1-1+pickidx)];
    PC0S0I1 = [PC0S0I1; PS0I1(end)*(1-sampleIdx_(var.C1-1+pickidx))];
    
%     <For monitoring>
%     if rem(nSample,1e3) == 0
%         meanS0 = sum(PS0)/nSample; varS0 = (1-meanS0)*meanS0/nSample; covS0 = sqrt(varS0)/meanS0;
%         meanC0S0I1 = sum(PC0S0I1)/sum(PS0I1); varC0S0I1 = (1-meanC0S0I1)*meanC0S0I1/sum(PS0I1);
%         covC0S0I1 = sqrt(varC0S0I1)/meanC0S0I1;
%         disp( [num2str( nSample ) '-th iter: E(S0) = ' num2str(meanS0) ' cov = ' num2str(covS0) ' | E(C0|S0I1) = ' ...
%             num2str(meanC0S0I1) ' cov = ' num2str(covC0S0I1)] )
%         eval( ['save demoSiouxFalls_MCS_' num2str(nSample/1e3)] );
%     end
end
    

%% Probability
meanS0 = sum(PS0)/nSample; varS0 = (1-meanS0)*meanS0/nSample; covS0 = sqrt(varS0)/meanS0;
meanC0S0I1 = sum(PC0S0I1)/sum(PS0I1); varC0S0I1 = (1-meanC0S0I1)*meanC0S0I1/sum(PS0I1);
covC0S0I1 = sqrt(varC0S0I1)/meanC0S0I1;
disp( [num2str( nSample ) '-th iter: E(S0) = ' num2str(meanS0) ' cov = ' num2str(covS0) ' | E(C0|S0I1) = ' ...
    num2str(meanC0S0I1) ' cov = ' num2str(covC0S0I1)] )

save demoSiouxFalls_MCS
