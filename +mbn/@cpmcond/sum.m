function newcpm = sum( cpm,vars,sumover )
%{
newcpm = SUM( cpm,vars,<sumover> )
sum cpm expelling vars or leaving vars

Input:
cpm: a cpm (cond)
vars: vars to be summed out or left
sumover: if 0 vars are summed out (default); 1 left

Output:
newcpm: resulted cpm

e.g.
cpm = cpmcond( 1:2,3:4,[1 1 2 3;1 2 2 3;2 1 2 2;-1 1 3 -1],[0.2 0.3 0.4 0.5] );
newcpm = sum( cpm,2 );
%}

import mbn.*

if nargin < 3
    sumover = 0;
end
if ~isempty( cpm.p )
    if sumover
        varout = setdiff( cpm.scope,vars );
    else
        if ~isempty( intersect(vars,cpm.scopep) )
            error( 'Conditioning variables cannot be summed out' )
        else
            varout = vars;
        end
    end

    [varleft,varleft_idx] = setdiff( cpm.scope,varout,'stable' );
    varleft = [varleft cpm.scopep];
    varleft_idx = [varleft_idx(:)' length(cpm.scope)+(1:length(cpm.scopep))];

    C = cpm.C; pold = cpm.p;
    Cold = C(:,varleft_idx);
    Cnew = []; pnew = [];
    while ~isempty(Cold)
        c_ = Cold(1,:);
        check = ismember( Cold,c_,'rows' );
        Cnew = [Cnew; c_];
        pnew = [pnew; sum(pold(check))];
        Cold(check,:) = []; pold(check) = [];
    end

    newcpm = cpm;
    newcpm.scope = setdiff( varleft,cpm.scopep,'stable' );
    newcpm.C = Cnew;
    newcpm.p = pnew;
else % if CPM is empty
    newcpm = cpm;
    switch sumover
        case 0
            newcpm.scope = setdiff( cpm.scope,vars,'stable' );
            newcpm.scopep = setdiff( cpm.scopep,vars,'stable' );
        case 1
            newcpm.scope = intersect( cpm.scope,vars,'stable' );
            newcpm.scopep = intersect( cpm.scopep,vars,'stable' );
    end
end