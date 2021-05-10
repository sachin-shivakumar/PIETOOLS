pvar s theta;

% find highest derivative in the PDE states
for i=0:max(der_pde)
        n_pde(i+1) = sum(der_pde==i);
end

PDE.n.no = no; PDE.n.nw = nw; PDE.n.no = nu; PDE.n.np = n_pde;
PDE.n.nz = nz; PDE.n.ny = ny; PDE.a = 0; PDE.b = 1;
PDE.var1 = s; PDE.var2 = theta;


N = length(n_pde)-1;

PDE.dyn.A = zeros(no);
PDE.dyn.Bxw = zeros(no,nw);
PDE.dyn.Bxu = zeros(no,nu);
PDE.dyn.Bxb = zeros(no,2*sum((0:N).*n_pde));


PDE.out.C1 = zeros(nz,no);
PDE.out.D1w = zeros(nz,nw);
PDE.out.D1u = zeros(nz,nu);
PDE.out.D1b = zeros(nz,2*sum((0:N).*n_pde));

PDE.out.C2 = zeros(ny,no);
PDE.out.D2w = zeros(ny,nw);
PDE.out.D2u = zeros(ny,nu);
PDE.out.D2b = zeros(ny,2*sum((0:N).*n_pde));

% PDE.out.Cv = [eye(no);zeros(nw,no);zeros(nu,no)];
% PDE.out.Dvw = [zeros(no,nw);eye(nw);zeros(nu,nw)];
% PDE.out.Dvu = [zeros(no,nu);zeros(nw,nu);eye(nu)];


np = sum(n_pde(1:N+1));
for i=1:N+1
    PDE.dyn.A0{i} = polynomial(zeros(np,sum(n_pde(i:N+1))));
    PDE.dyn.A1{i} = polynomial(zeros(np,sum(n_pde(i:N+1))));
    PDE.dyn.A2{i} = polynomial(zeros(np,sum(n_pde(i:N+1))));
    PDE.dyn.Bxr{i} = polynomial(zeros(no,sum(n_pde(i:N+1))));
    PDE.out.Dzr{i} = polynomial(zeros(nz,sum(n_pde(i:N+1))));
    PDE.out.Dyr{i} = polynomial(zeros(ny,sum(n_pde(i:N+1))));
    PDE.BC.Br{i} = polynomial(zeros(sum((0:N).*n_pde),sum(n_pde(i:N+1))));
end

PDE.dyn.Dv = polynomial(zeros(np,no+nw+nu));
PDE.dyn.Db = polynomial(zeros(np,2*sum((0:N).*n_pde)));

PDE.BC.B = zeros(sum((0:N).*n_pde), 2*sum((0:N).*n_pde));
PDE.BC.Bx = zeros(sum((0:N).*n_pde),no);
PDE.BC.Bw = zeros(sum((0:N).*n_pde),nw);
PDE.BC.Bu = zeros(sum((0:N).*n_pde),nu);