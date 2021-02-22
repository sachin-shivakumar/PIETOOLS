dim = Top.dim;

%get dim of input vectors
m = dim(1,2); n = dim(2,2);
N=5;
%
x_f = polynomial(zeros(m+n,1)); 

%generate pvar coeffs
if m>0
for i=1:m
    eval(['pvar cf_' num2str(i)]);
    x_f(i,1) = eval(['cf_' num2str(i)]);
end
for i=1:n
    for j =0:N
        eval(['pvar cf_' num2str(m+1+j+(i-1)*N)]);
        x_f(m+i,1) = x_f(m+i,1) + eval(['cf_' num2str(m+1+j+(i-1)*N)])*cheb_poly_2(s,j);
    end
end
elseif m==0
for i=1:n
    for j =0:N
        eval(['pvar cf_' num2str(j+(i-1)*N)]);
        x_f(i,1) = x_f(i,1) + eval(['cf_' num2str(j+(i-1)*N)])*cheb_poly_2(s,j);
    end
end
end

nw = B1op.dim(1,2); nu = B2op.dim(1,2);
w_f = eye(nw); u_f = eye(nu);

Topd = discretize_opvar(Top,N,x_f);
Aopd = discretize_opvar(Aop,N,x_f);
% B1opd = discretize_opvar(B1op,N,w_f);
% B2opd = discretize_opvar(B2op,N,u_f);

Tbar = Topd.coefficient; Abar = Aopd.coefficient;
c0 = rand(6,1);
dydt = @(t,y) Abar*y;

options = odeset('Mass',Tbar);
[t,y] = ode15s(dydt, [0,20], c0,options);

dx = linspace(-1,1,50);
x_sol = double(subs(subs(x_f,Topd.varname,y'),s,dx));
