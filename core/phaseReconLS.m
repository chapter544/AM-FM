function y = phaseReconLS(U, V, Pinit, PHASE_CONST) 

U = PHASE_CONST .* U;
V = PHASE_CONST .* V;

P = poisson_solver_function_neumann(U, V);

% Equalize the biased constant
y = P + ones(size(P)) .* (Pinit - P(1,1));
