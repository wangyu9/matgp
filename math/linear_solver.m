function [x] = linear_solver(A,B,method)

switch(method)
    case 'cholmod'
        [x,stats] = cholmod2(A,B);
        fprintf('SuiteSparse Memory: %f.\n',stats(5));
    case 'default'
        x = A\B;
    otherwise
        fprintf('Unknown linear solver type, use default instead.\n');
        x = A\B;
end