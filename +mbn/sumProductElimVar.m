function newCpm = sumProductElimVar( Cpm,elimvar )
%{
newCpm = SUMPRODUCTELIMVAR( Cpm,elimvar )
Eliminate a rv from a set of CPMs by sum-product

Input:
Cpm: a set of CPMs (cell array)
elimvar: rv to be eliminated
--
Output:
newCpm: new set of CPMs w/o var

e.g. (Student BN in PGM (2010))
d=1;i=2;g=3;
Cpm{d} = cpmjoint( d,[0; 1],[0.6 0.4] );
Cpm{i} = cpmjoint( i,[0; 1],[0.7 0.3] );
Cpm{g} = cpmcond( g,[i d],[repmat((1:3)',4,1) repelem((0:1)',6,1) repmat(repelem((0:1)',3,1),2,1)], ...
 [0.3 0.4 0.3 0.05 0.25 0.7 0.9 0.08 0.02 0.5 0.3 0.2]);
elimvar = d;
newCpm = sumProductElimVar( Cpm,elimvar )
%}

import mbn.*
check = zeros( size(Cpm) );
for ii = 1:length(Cpm)
    cpm_ = Cpm{ii};
    if isa( cpm_,'cpmcond' )
        check(ii) = ismember( elimvar,[cpm_.scope cpm_.scopep] );
    else
        check(ii) = ismember( elimvar,cpm_.scope );
    end
end

check = logical(check);
Cpm1 = Cpm(check); newCpm = Cpm(~check);

if ~isempty(Cpm1)
    cpm_ = Cpm1{1};
    for ii = 2:sum(check)
        cpm_ = product(cpm_,Cpm1{ii});
    end
    cpm_ = sum(cpm_,elimvar);
    newCpm = [newCpm(:)' {cpm_}];
end
