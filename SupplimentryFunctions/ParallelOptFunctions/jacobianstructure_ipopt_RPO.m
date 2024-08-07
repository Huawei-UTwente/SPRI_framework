function Js = jacobianstructure_ipopt_RPO(auxdata)
% sparse matrix of the jacobian

Js = sparse(auxdata.row, auxdata.col, ones(1, length(auxdata.row)));

end