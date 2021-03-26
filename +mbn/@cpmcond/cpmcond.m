classdef cpmcond
% Conditional pmf: P(scope|scopep)
    properties
        scope;
        scopep;
        C;
        p;
        val;
    end
    
    methods
        function cpmc = cpmcond(scope, scopep, C, p, val)
            % cpmcond(<scope>, <scopep>, <C>, <p>, <val>)
            % eg p1 = cpmcond(3,[1 2],[1 1 0;1 2 1;1 2 2;2 0 2],[0.5 0.3 0.2 1])
            % or p1=cpmcond; p1.scope = 3; p1.scopep=[1 2]; p1.C=[1 1 0;1 2 1;1 2 2;2 0 2]; p1.p=[0.5 0.3 0.2 1];
            if nargin>0
                if ~isnumeric(scope)
                    error( 'Scope1 must be a numerical vector' )
                else
                    cpmc.scope = scope(:)';
                end
                
                if nargin>1
                    if ~isnumeric(scopep)
                        error( 'Scope2 must be a numerical vector' )
                    elseif ~isempty(intersect(scope,scopep))
                        error('Scope 1 and 2 must not have common vars')
                    else
                        cpmc.scopep = scopep(:)';
                    end
                    
                    if nargin>2
                        if ~isnumeric(C)
                            error('Event matrix must be a numerical array')
                    elseif size(C,2) ~= length([scope(:)' scopep(:)'])
                        error( 'The # of cols in C must be same with the # of vars in scope' )
                    else
                        cpmc.C = C;
                        end
                    
                        if nargin>3
                            if ~isnumeric(p)
                                error( 'Prob vector must be a numerical vector' )
                            elseif ~isvector(p) && ~isempty(p)
                                error('Prob vector must be a vector')
                            elseif size(C,1)~=length(p) && ~isempty(p) && ~isempty(C)
                                error('The # of rules in C must be same with that in p')
                            else
                                cpmc.p = p(:); 
                            end
                            
                            if nargin>4
                                cpmc.val = val;
                            end
                        end
                    end
                end
            end
        end
    end
end