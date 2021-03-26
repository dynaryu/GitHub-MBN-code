% Construct BN for Sioux Falls Network

clear; close all;
import mbn.*
import Code_Sioux_fragility.*
%% Info.
SFnet.locnode = load('data/SiouxFalls_node.txt');
SFnet.locnode = SFnet.locnode(:,2:3)*0.000025; % unit :in -> km 
SFnet.arcinfo = load('data/SiouxFalls_net1.txt');
SFnet.locarc = ( SFnet.locnode( SFnet.arcinfo(:,1),: ) + SFnet.locnode( SFnet.arcinfo(:,2),: ) ) * .5;
EQloc = 0.5*[(-4:5)' (5:-1:-4)']; % EQ location (x,y)
% EQloc = 0.5*[(-6:3)' (3:-1:-6)']; % EQ location (x,y)

run( 'SiouxFalls_fig.m' )

%% Variables
nComp = size(SFnet.arcinfo,1);

% Magnitude: X1
var.M = 1; 
Mvars = linspace(6,8.5,6);
cpmSioux{var.M} = cpmjoint;
cpmSioux{var.M}.scope = var.M;
cpmSioux{var.M}.C = ( 1:(length(Mvars)-1) )';
cpmSioux{var.M}.p = 1 / ( 1-exp( -0.76*(8.5-6) ) ) * (1-exp(-0.76*(Mvars-6)));
cpmSioux{var.M}.p = diff(cpmSioux{var.M}.p);
cpmSioux{var.M}.p = cpmSioux{var.M}.p(:);
cpmSioux{var.M}.val = (Mvars(1:end-1) + Mvars(2:end) )' *.5;

% EQ location: X2
var.L = 2; 
cpmSioux{var.L} = cpmjoint;
cpmSioux{var.L}.scope = var.L;
cpmSioux{var.L}.C = (1:10)';
cpmSioux{var.L}.p = ones(length(cpmSioux{var.L}.C),1) / length(cpmSioux{var.L}.C);
cpmSioux{var.L}.val = EQloc; % (x,y)

% Deterioration: X3~X78
ke.coef = [0.924 0.265 0.676];
[M,V] = lognstat(1+log(1e2./ke.coef),0.05);
ke.prob = normcdf(log(20),M,V);

scen{1,1} = [39 37 38 74 66 75 73 76 62 64 65 69 70 72 42 71]; % Tidal: ke = 0.924
scen{2,1} = [59 61 63 68 46 67 57 45 44 41 34 40 33 36 7 35]; % Splash: ke = 0.265
scen{3,1} = setdiff( 1:nComp,union(scen{1},scen{2}) ); % atmospheric: ke = 0.676

for ii = 1:nComp
    eval(['var.D' num2str(ii) '=' num2str(ii+var.L) ';']);
    tmp_ = find( cellfun( @(x) ismember(ii,x),scen ) );
    cpmSioux{ii+var.L} = cpmjoint( ii+var.L,[0 1]',[1-ke.prob(tmp_) ke.prob(tmp_)]',{'False' 'True'}' ); % 1: corr O, 0: corr X
end

% Inspection: X79~158
eval( ['d_ = var.D' num2str(nComp) ';'] );
for ii = 1:nComp
    eval(['var.I' num2str(ii) '=' num2str(ii+d_) ';']);
    cpmSioux{ii+d_} = cpmcond( ii+d_,d_-nComp+ii,[0 0;1 0;0 1;1 1],[0.9 0.1 0.2 0.8],{'False' 'True'}' ); % 1: corr detected, 0: Not detected
end

% Component Capacity
load SiouxFalls_BN % Lee et al. (2011)'s example code has been used for bridges' reliabilities

% System
load SiouxFalls_RDA pairs
SF_RDA = load('SiouxFalls_RDA','G');
SF_RDA = SF_RDA.G;
cq_ = {'Fail','Survive'}; % 'Connected': IN,MD | 'Disconnected': HV,CP
c_ = 230;
var.S = 231;
nGraph = length(SF_RDA.PATHNODE); % graphs to be considered 
sc_ = -ones(nGraph,nComp+1); sp_ = ones(nGraph,1);
for ii = 1:nGraph
    sur_ = union(SF_RDA.PATHLINK{ii}, SF_RDA.S{ii});
    fail_ = SF_RDA.S_{ii};
    if ~isempty( intersect(sur_,fail_) ); error( 'Conflicting Components!' ); end
    sc_(ii,sur_+1) = 1;
    sc_(ii,fail_+1) = 0;    
    
    if ~isempty(SF_RDA.PATHLINK{ii})
        sc_(ii,1) = 1;
    else
        sc_(ii,1) = 0;
    end
end
cpmSioux{var.S} = cpmcond( var.S,c_-nComp+1:c_,sc_,sp_,cq_ );

save SiouxFalls_BN.mat Graph nComp cpmSioux SFnet var pairs