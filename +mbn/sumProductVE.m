function newcpm = sumProductVE( Cpm,elimvars )
%{
newCpm = SUMPRODUCTVE( Cpm,elimvars,order )
sum-product-variable-elimination

Input:
Cpm: a cell array of CPMs
elimvars: variables to be eliminated (in elimination ordering)
--
Output:
newcpm: a CPM after elimination

e.g. (Student BN in PGM (2010))
d=1;i=2;g=3;s=4;l=5;a=6;j=7;
Cpm{d} = cpmjoint( d,[0; 1],[0.6 0.4] );
Cpm{i} = cpmjoint( i,[0; 1],[0.7 0.3] );
Cpm{g} = cpmcond( g,[i d],[repmat((1:3)',4,1) repelem((0:1)',6,1) repmat(repelem((0:1)',3,1),2,1)], ...
 [0.3 0.4 0.3 0.05 0.25 0.7 0.9 0.08 0.02 0.5 0.3 0.2]);
Cpm{s} = cpmcond(s,i,[0 0;1 0;0 1;1 1],[0.95 0.05 0.2 0.8]);
Cpm{l} = cpmcond(l,g,[repmat((0:1)',3,1) repelem((1:3)',2,1)],[0.1 0.9 .4 .6 .99 .01]);
Cpm{a} = cpmjoint(a,[0;1],[.3 .7]);
Cpm{j} = cpmcond(j, [a s l], [0 0 -1 -1;1 0 -1 -1;0 1 1 -1;1 1 1 -1;0 1 0 0;1 1 0 0;0 1 0 1;1 1 0 1], ...
[0.8 0.2 0.1 0.9 .9 .1 .4 .6]);
elimvars = [d i s l];
newcpm = sumProductVE( Cpm,elimvars );
%}

import mbn.*

newCpm = Cpm;
for ii = 1:length(elimvars)
    newCpm = sumProductElimVar( newCpm,elimvars(ii) );
end
newcpm = newCpm{1};
for ii = 2:length(newCpm)
    newcpm = product( newcpm,newCpm{ii} );
end