function [compatcheck,reducedC] = compat( scopec,c,scopeC,C,assignsetting )
%{
COMPATCHECK = COMPAT( scopec,c,scopeC,C )
check the compatibility b/w an assignment c and a set of assignments C

Input:
scopec: scope of c (array)
c: an assignment c
scopeC: scope of C (array)
C: a set of assignments
assignsetting: if 0 the assignments in reducedC NOT changed (default); 1 "-1" state changed to the values in scopec
---
Output:
compatcheck: 1 x size(C,1) array with 1 if compatible w/ c; 0 otherwise
reducedC: C containing only the rules compatible with c

e.g.
compatcheck = compat( [2 1 3],[1 3 -1],[3 1 4],[1 3 1;2 3 2;3 -1 3;-1 3 2;-1 2 1;-1 1 1] );
%}

import mbn.*
% if ~isvector(scopec) || ~isvector(scopeC)
%     error('Scope must be a numerical vector')
% end
% if ~isvector(c)
%     error('1st input rule must be a vector')
% end

idxnegc = find( c==-1 );
[~,idxc] = intersect(scopec,scopeC,'stable');
idxc = setdiff(idxc,idxnegc,'stable');
[~,idxC] = ismember(scopec(idxc),scopeC);

compatcheck_ = 1:size(C,1);
if nargout<2
    for ii = 1:length(idxc)
        C_ = C(compatcheck_,:);
        tmp_ = union( find( C_(:,idxC(ii)) == c( idxc(ii) ) ), find( C_(:,idxC(ii))==-1 ) );
        compatcheck_ = compatcheck_(tmp_);
    end
    compatcheck = zeros(size(C,1),1); compatcheck(compatcheck_) = 1;
    compatcheck = logical(compatcheck);
elseif nargout<3
    if nargin<5
        assignsetting = 0;
    end
    compatcheck_ = 1:size(C,1);
    if assignsetting
        for ii = 1:length(idxc)
            C_ = C(compatcheck_,:);
            compatcheck1 = find( C_(:,idxC(ii)) == c( idxc(ii) ) );
            compatcheck2 = find( C_(:,idxC(ii))==-1 );
            C(compatcheck_(compatcheck2),idxC(ii)) = c( idxc(ii) );
            compatcheck_ = union( compatcheck_(compatcheck1),compatcheck_(compatcheck2) );
        end
    else
        for ii = 1:length(idxc)
            C_ = C(compatcheck_,:);
            tmp_ = union( find( C_(:,idxC(ii)) == c( idxc(ii) ) ), find( C_(:,idxC(ii))==-1 ) );
            compatcheck_ = compatcheck_(tmp_);
        end
    end
    compatcheck = zeros(size(C,1),1); compatcheck(compatcheck_) = 1;
    compatcheck = logical(compatcheck);
    reducedC = C( compatcheck_,: );
end