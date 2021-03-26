%{
Get the approximate inference result of P(C1=1|S=1)
in the example RBD
%}

clear; close all;
import mbn.*

nComp = 8; % # of C's
var.s = nComp+1;
%% [1]determinstic bounds

% problem
for ii = 1:nComp
    cpm{ii,1} = cpmjoint( ii,[0;1],[0.1;0.9],{'Fail';'Survive'} );
end
cpm{var.s} = cpmcond( var.s,1:nComp );
C_ = -ones(7,nComp+1);
C_(1,[1 9]) = [0 0];
C_(2,[1 8 9]) = [0 0 1];
C_(3,[1 2 8 9]) = [1 1 1 1];
C_(4,[1 2 3 8 9]) = [1 0 1 1 1];
C_(5,[1 2 3 4 8 9]) = [1 0 0 1 1 1];
C_(6,[1 2 3 4 5 8 9]) = [0 0 0 0 0 1 1];
C_(7,[1 2 3 4 5 6 8 9]) = [0 0 0 0 0 0 1 1];
cpm{var.s}.C = C_; cpm{var.s}.p = ones(size(C_,1),1); 

% Inference
elimvars = 2:nComp;
cpmSC1 = sumProductVE(cpm,elimvars);

% P(S=1,C1=1)
c_ = ismember( cpmSC1.C,[1 1],'rows' );
Ps1c1_lower = sum( cpmSC1.p(c_) );
Ps1c1_upper = 1-sum( cpmSC1.p(~c_) );
disp( ['The bounds of P(S=1,C1=1) is [ ',num2str(Ps1c1_lower), ' ', num2str(Ps1c1_upper),' ]'] )

% P(S=1)
cpmS = sum(cpmSC1,var.s,1);
c_ = ismember( cpmS.C,1,'rows' );
Ps1_lower = sum( cpmS.p(c_) );
Ps1_upper = 1-sum( cpmS.p(~c_) );
disp( ['The bounds of P(S=1) is [ ',num2str(Ps1_lower), ' ', num2str(Ps1_upper),' ]'] )

% P(C1=1|S=1)
disp( ['The bounds of P(C1=1|S=1) is [ ',num2str(Ps1c1_lower/Ps1_upper), ' ', num2str(Ps1c1_upper/Ps1_lower),' ]'] )

%% [2]sampling
rng(77)
nSample = 10;
PC1 = .9;
Csample = rand(nSample,nComp)<PC1;
Q = PC1.^sum(Csample,2) .* (1-PC1).^( nSample-sum(Csample,2) );
Ssample = zeros(nSample,1);
for ii = 1:nSample
    switch Csample(ii,8); case 0; Ssample(ii) = 0;
        case 1
            switch Csample(ii,7); case 0; Ssample(ii) = 0;
                case 1
                    switch Csample(ii,1); case 1; Ssample(ii) = 1;
                        case 0
                            switch Csample(ii,2); case 1; Ssample(ii) = 1;
                                case 0
                                    switch Csample(ii,3); case 1; Ssample(ii) = 1;
                                        case 0
                                            switch Csample(ii,4); case 0; Ssample(ii) = 0;
                                                case 1
                                                    switch Csample(ii,5); case 0; Ssample(ii) = 0;
                                                        case 1
                                                            switch Csample(ii,6); case 0; Ssample(ii) = 0;
                                                                case 1; Ssample(ii) = 1;
                                                            end
                                                    end
                                            end
                                    end
                            end
                    end
            end
    end
end
D.w = Ssample; % P=Q
D.mud = Ssample.* Csample(:,1);
PC1S1 = sum( D.w.*D.mud ) / sum(D.w);
PC1S1_var = sum( D.w.^2 .*( D.mud-PC1S1 ).^2 ) / sum(D.w)^2;
PC1S1_cov = sqrt(PC1S1_var/sum(Ssample))/PC1S1;
disp( ['The bounds of P(C1=1|S=1) is [ ',num2str(Ps1c1_lower/Ps1_upper), ' ', num2str(Ps1c1_upper/Ps1_lower),' ]'] )

