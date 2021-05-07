function PDE_text = executive_pde_text_cleaner(Text, PDE_vars, ODE_vars, Dis_vars, Inp_vars, RegOut_vars, ObsOut_vars)

PDE_text = erase(Text,'$'); % can be dynamics or BCs
PDE_text = strrep(PDE_text, '\theta', 'Theta');


for i=1:size(PDE_vars,1)
    PDE_text = strrep(PDE_text, PDE_vars(i,1), ['X' num2str(i)]);
end

for i=1:size(ODE_vars,1)
    PDE_text = strrep(PDE_text, ODE_vars(i,1), ['x' num2str(i)]);
end

for i=1:size(Dis_vars,1)
    PDE_text = strrep(PDE_text, Dis_vars(i,1), ['w' num2str(i)]);
end

for i=1:size(Inp_vars,1)
    PDE_text = strrep(PDE_text, Inp_vars(i,1), ['u' num2str(i)]);
end

for i=1:size(RegOut_vars,1)
    PDE_text = strrep(PDE_text, RegOut_vars(i,1), ['z' num2str(i)]);
end

for i=1:size(ObsOut_vars,1)
    PDE_text = strrep(PDE_text, ObsOut_vars(i,1), ['y' num2str(i)]);
end
end