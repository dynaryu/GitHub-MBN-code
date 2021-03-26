function newcpm = sum( cpm,vars,sumover )
%{
newcpm = SUM( cpm,vars,<sumover> )
sum cpm expelling vars or leaving vars

Input:
cpm: a cpm (joint)
vars: vars to be summed out or left
sumover: if 0 vars are summed out (default); 1 left

Output:
newcpm: resulted cpm

e.g.
cpm = cpmjoint( 1:4,[1 1 2 3;1 2 2 3;2 1 2 2;-1 1 3 -1],[0.2 0.3 0.4 0.5] );
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
        varout = vars;
    end

    [varleft,varleft_idx] = setdiff( cpm.scope,varout,'stable' );

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
    newcpm.scope = varleft;
    newcpm.C = Cnew;
    newcpm.p = pnew;
else % if CPM is empty
    newcpm = cpm;
    switch sumover
        case 0
            newcpm.scope = setdiff( cpm.scope,vars,'stable' );
        case 1
            newcpm.scope = intersect( cpm.scope,vars,'stable' );
    end
end