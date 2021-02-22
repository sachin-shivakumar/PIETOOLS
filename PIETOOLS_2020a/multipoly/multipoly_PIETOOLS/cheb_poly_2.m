function Un = cheb_poly_2(var, N)
% This function returns a chebyshev polynomial of the Nth degree in the
% polynomial variables var

U0 =1;
U1 = var;

Un2 = U0;
Un1 = U1;

if N==0
    Un = U0;
elseif N==1
    Un = U1;
else
    for i=2:N
        Un = 2*var*Un1 - Un2;
        Un2 = Un1;
        Un1 = Un;
    end
end
end