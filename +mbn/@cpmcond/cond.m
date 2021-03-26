function newcpm = cond( cpm,E,ev,probsetting )
%{
newcpm = CONDITIONING( cpm,E,ev,<probsetting> )

Input:
cpm: a cpm (cond)
E: set of ev
ev: given context w/o "-1" state
probsetting: 0 - p NOT changed (default), 1: p changed to 1 (in case ev==c)
--
Output:
newcpm: contains only the rules compatible w/ E=ev, with original or changed p

e.g.
cpm = cpmcond( 1:2,3:4,[1 1 2 3;1 2 -1 3;2 1 2 2;-1 1 3 -1],[0.2 0.3 0.4 0.5] ); E = [1 2 4]; ev = [1 1 3];
newcpm = conditioning( cpm,E,ev,1 );
%}

import mbn.*
if sum( ev==-1 )
    error( 'Given context must not contain "-1" state' )
end

scope = [cpm.scope cpm.scopep];

[scope_,idxc] = intersect(scope,E,'stable');
[~,idxe] = ismember(scope_,E);

C_ = cpm.C(:,idxc);
ev_ = ev( idxe );
[compatcheck,Cnew_] = compat( scope_,ev_,scope_,C_,1 );
Cnew = cpm.C(compatcheck,:);
Cnew( :,idxc ) = Cnew_;

if nargin>3
    if probsetting
        if isempty( setdiff(cpm.scope,E) )
            pnew = ones( sum(compatcheck),1 );
        end
    else
        pnew = cpm.p(compatcheck);
    end
else
    pnew = cpm.p(compatcheck);
end

newcpm = cpm; newcpm.C = Cnew; newcpm.p = pnew;