%build the mass and stiffness matrices that describe the 2nd order system.
%INPUTS
%string_params: a struct containing the system parameters describing the string
% string_params.n: number of masses
% string_params.M: total mass attached to the string
% string_params.Uf_func: function describing motion of end point
% string_params.dUfdt_func: time derivative of Uf
% string_params.Tf: %tension in string
% string_params.L: %length of string
% string_params.c: %damping coefficient
% string_params.dx: %horizontal spacing between masses
%OUTPUTS
%M_mat: the n x n mass (inertia) matrix
%K_mat: the n x n stiffness matrix
function [M_mat,K_mat] = construct_2nd_order_matrices(string_params)
    M_total=string_params.M;
    n=string_params.n;
    Tf=string_params.Tf;
    dx=string_params.dx;
    I=eye(n);
    Q=-2*I+circshift(I,[0,1])+circshift(I,[0,-1]);
    Q(1,end)=0;  Q(end,1)=0;

    M_mat=M_total/n*I;
    K_mat=-Tf/dx*Q;
end
