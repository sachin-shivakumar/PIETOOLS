function [L] = getObserver(P,Z)
% This function returns the controller gains K for a PIE system given lyapunov
% operator P and PI operator Z

if isvalid(P) && isvalid(Z)
L = inv_opvar(P)*Z;
else
    error("Inputs must be opvar variables");
end
end