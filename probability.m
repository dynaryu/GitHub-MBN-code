function prob = probability( scopec,c,cpm )
%{
prob = PROBABILITY( scopec,c,cpm )
get the probability of context c with scopec

INPUT:
scopec: scope of given context (when matrices --> each scope for each context; row vector --> all given contexts share scope)
c: rows of given context
cpm: a cpm for which prob is evaluated.
--
OUTPUT:
prob: probability of c w.r.t. cpm
%}
import mbn.*
%%
prob = zeros( size(c,1),1 );
if size(scopec,1) == 1
    for ii = 1:length(prob)
        tmp = compat( scopec, c(ii,:), cpm.scope, cpm.C );
        prob(ii) = sum( cpm.p(tmp) );
    end
else
    for ii = 1:length(prob)
        tmp = compat( scopec(ii,:), c(ii,:), cpm.scope, cpm.C );
        prob(ii) = sum( cpm.p(tmp) );
    end
end
