function newcpm = product( cpm1,cpm2 )
%{
newcpm = PRODUCT( cpm1,cpm2 )
Product cpm1 and cpm2

Input:
cpm1, cpm2: cpm

Output:
newcpm: cpm1*cpm2

e.g.
cpm1=cpmjoint([1 3],[1 1;2 1;1 2;2 2],[0.1 0.9 0.4 0.6]); cpm2=cpmjoint([2 4],[1 1;2 1;1 2;2 2],[0.3 0.7 0.9 0.1]);
or
cpm1 = cpmcond( 1,[2 3],[1 1 -1;2 1 -1;1 2 1;2 2 1],...
[0.2 0.8 0.3 0.7] );
cpm2 = cpmcond( [3 4],1,[1 1 1;2 1 1;1 1 2;2 1 2],[0.8 0.2 0.3 0.7] );
%}

import mbn.*
isj = [isa(cpm1,'cpmjoint') isa(cpm2,'cpmjoint')];

% error check
scope1_ = cpm1.scope; scope2_ = cpm2.scope;
if ~isempty( intersect(scope1_,scope2_) )
    error( 'CPMs must not have common variables to be multiplied' )
end

% 
if ~isj(1); scope1_ = [cpm1.scope cpm1.scopep]; end
if ~isj(2); scope2_ = [cpm2.scope cpm2.scopep]; end
    
comvars = intersect( scope1_,scope2_,'stable');
[~,com_idx1] = ismember(comvars,scope1_);
[~,com_idx2] = ismember(comvars,scope2_);
[~,marg_idx1] = setdiff(scope1_,comvars,'stable');
[~,marg_idx2] = setdiff(scope2_,comvars,'stable');
scope_ = unique( [comvars scope1_ scope2_],'stable' );

C1_com = cpm1.C(:,com_idx1);
C2_com = cpm2.C(:,com_idx2);
Cnew = []; pnew = [];
for ii = 1:size(C1_com,1)
    mu1_ = C1_com(ii,:);
    posidx1_ = ( mu1_~=-1 );
    [check_,Cnew_] = compat( comvars(posidx1_),mu1_(posidx1_),comvars,C2_com,1 );
    Cnew_ = [Cnew_ repmat(cpm1.C(ii,marg_idx1),sum(check_),1) cpm2.C(check_,marg_idx2)];
    pnew_ = log( cpm1.p(ii) ) + log( cpm2.p(check_) );
    pnew_ = exp(pnew_);
    Cnew = [Cnew; Cnew_]; pnew = [pnew; pnew_];
end

% 
if sum(~isj)<2
    scope = comvars; scopep = [];
else % both are cpmcond
    scopep = intersect(cpm1.scopep, cpm2.scopep,'stable');
    scope = setdiff( comvars,scopep,'stable' );
end
if ~isj(1) 
    scope = [scope setdiff(cpm1.scope,comvars,'stable')];
    scopep = [scopep setdiff(cpm1.scopep,comvars,'stable')];
    if ~isj(2)
        scope = [scope setdiff(cpm2.scope,comvars,'stable')];
        scopep = [scopep setdiff(cpm2.scopep,comvars,'stable')];
    else
        scope = [scope setdiff(cpm2.scope,comvars,'stable')];
    end
else
    scope = [scope setdiff(cpm1.scope,comvars,'stable')];
    if ~isj(2)
        scope = [scope setdiff(cpm2.scope,comvars,'stable')];
        scopep = [scopep setdiff(cpm2.scopep,comvars,'stable')];
    else
        scope = [scope setdiff(cpm2.scope,comvars,'stable')];
    end
end

if isempty(cpm1.p) + isempty(cpm2.p)
    Cnew = zeros(0,length([scope scopep])); pnew = zeros(0,1);
else    
    [~,idx_] = ismember( [scope scopep],scope_ );
    Cnew = Cnew(:,idx_);
end

if ~isempty( scopep )        
    newcpm = cpmcond( scope,scopep,Cnew,pnew );
else
    newcpm = cpmjoint( scope,Cnew,pnew );
end
    