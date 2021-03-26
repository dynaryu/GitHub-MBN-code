% Get a subset of rules from RDA in Sioux Falls network example

clear; close all;
import mbn.*

%%
pairs = load('data/SiouxFalls_net1.txt');
nNode = max(pairs(:));
nLink = size(pairs,1);
Adj = zeros(nNode); Anumber = zeros(nNode);
for ii = 1:nLink
    Adj(pairs(ii,1),pairs(ii,2)) = 1;
    Anumber(pairs(ii,1),pairs(ii,2)) = ii;
end
Prob = Adj;% dummy
Prob(Prob>0) = 0.1; 
ProbLink = Prob(Prob>0);% dummy (leading to original RDA, not sRDA)
sNode = 13; tNode = 2;

[G,sysR,sysF,err] = sRDA(Adj,Anumber,Prob,ProbLink,sNode,tNode,0,1,5e3);

save SiouxFalls_RDA Adj Anumber G nLink nNode pairs Prob ProbLink sNode tNode sysF sysR
%%