function newcpm = conditioning( cpm,E,ev,probsetting )
%{
newcpm = CONDITIONING( cpm,E,ev,<probsetting> )

Input:
cpm: either a cpm or a cell array of cpm
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

if nargin<4
    probsetting = 0;
end
if ~iscell( cpm )
    newcpm = cond( cpm,E,ev,probsetting );
else
    newcpm = cell( size(cpm) );
    for ii=1:length(cpm)
        newcpm{ii} = cond( cpm{ii},E,ev,probsetting );
    end
end