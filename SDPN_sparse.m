function [ D ] = SDPN_sparse( A, N, E, beta, gamma,  source )
% INPUT: adjacency matrix A (unweighted), number of nodes N, number of edges E,
% parameters of SIR beta,gamma, id of source node
% OUTPUT: Vector of first infection times from source to other nodes

Rho = exprnd(1/beta, E, 1);

Tau_node = exprnd(1/gamma,N,1);
[I,J]=find(A);
Tau = Tau_node(I);

W = sparse(I,J,Rho.*(Rho <= Tau),N,N);
D = graphshortestpath(W, source);

end

