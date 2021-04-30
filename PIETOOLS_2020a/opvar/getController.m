function [K] = getController(P,Z)
% This function returns the controller gains K for a PIE system given lyapunov
% operator P and PI operator Z

if isvalid(P) && isvalid(Z)
K = Z*inv_opvar(P);
else
    error("Inputs must be opvar variables");
end
end