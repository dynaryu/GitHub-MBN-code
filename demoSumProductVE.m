%{
Demonstrate the use of CPM in Sum-Product-VE algorithm
Example: Student BN in PGM (2010)
%}

close all; clear
import mbn.*
%% Problem
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

%% Inference
elimvars = [d i s l];
newcpm = sumProductVE( Cpm,elimvars );

j1g2a1 = compat( [j g a],[1 2 1],newcpm.scope,newcpm.C );
g2a1 = compat( [g a],[2 1],newcpm.scope,newcpm.C );
disp( ['P(j1|g2,a1) = ' num2str(sum(newcpm.p(j1g2a1))/sum(newcpm.p(g2a1)))] )