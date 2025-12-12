%modal analysis

function string_simulation_template03()
    num_masses = 200; %your code here
    total_mass = 10; %your code here
    tension_force = 5; %your code here
    string_length = 20; %your code here
    damping_coeff = 0.0067; %your code here
    dx = string_length/(num_masses+1);

    string_params = struct();
    string_params.n = num_masses;
    string_params.M = total_mass;
    string_params.Tf = tension_force;
    string_params.L = string_length;
    string_params.c = damping_coeff;
    string_params.dx = dx;
    
    mode_index=4; 

    rho = total_mass/string_length;
    c=sqrt(tension_force/rho);
    w_pulse=string_length/20;
    h_pulse=7;
    

    %[M_mat,K_mat] = construct_2nd_order_matrices(string_params);
    % [Ur_mat,lambda_mat] = eig(K_mat,M_mat); 

    % mode_shape_LA = Ur_mat(:,mode_index);

    %omega_n=sqrt(lambda_mat(mode_index,mode_index));
    omega_n_spatial=pi*mode_index/string_length;
    omega_n=c*omega_n_spatial;

    xlist = linspace(0,string_length,num_masses+2);
    xlist = xlist(2:end-1);

    mode_shape_WE=sin(omega_n_spatial*xlist);

    amplitude_Uf = 3;
    omega_Uf = omega_n;

    %mode = Ur_mat(:,mode_index);

    Uf_func = @(t_in) amplitude_Uf*cos(omega_Uf*t_in);
    dUfdt_func = @(t_in) -omega_Uf*amplitude_Uf*sin(omega_Uf*t_in);

    string_params.Uf_func = Uf_func;
    string_params.dUfdt_func = dUfdt_func;

     

    %U0 = triangle_pulse(xlist,w_pulse, h_pulse)';   
    %dUdt0 = -c*triangle_pulse_derivative(xlist,w_pulse, h_pulse)';
    U0 = zeros(num_masses,1);     
    dUdt0 = zeros(num_masses,1);
    V0 = [U0;dUdt0];

    my_rate_func = @(t_in,V_in) string_rate_func01(t_in,V_in,string_params);

    tspan = linspace(0,3*string_length/c,5001);
    [t_list,V_list] = ode45(my_rate_func,tspan,V0);
    V_list = V_list';
    
    xlist=[0,xlist,string_length];
    mode_shape_WE=[0,mode_shape_WE,0];
    hold on;
    mode_shape_plot=plot(0,0,'o-','color','b','LineWidth',2,'markerfacecolor','k','markeredgecolor','k','markersize',5);
    set(mode_shape_plot,'xdata', xlist, 'ydata', mode_shape_WE)

    ball_plot_struct = initialize_balls_plot();
    axis([0,string_length,-5,5]);
    xlabel("x")
    ylabel("f(x)")

    for k = 1:length(t_list)
        update_balls_plot(ball_plot_struct,V_list(:,k),t_list(k),string_params);
        drawnow;
    end
end

%triangle pulse function
%INPUTS:
%t: current time
%w: width of pulse (starts at t=0, ends at t=h)
%h: height of pulse
%OUTPUTS:
%res: pulse evaluated at t
function res = triangle_pulse(t,w,h)
    t = t*(2/w);
    res = 1-min(1*abs(t-1),1);
    res = h*res;
end
%triangle pulse function (derivative)
%INPUTS:
%t: current time
%w: width of pulse (starts at t=0, ends at t=h)
%h: height of pulse
%OUTPUTS:
%res: derivative of pulse evaluated at t
function res = triangle_pulse_derivative(t,w,h)
    t = t*(2/w);
    res = -sign(t-1).*(abs(t-1)<1);
    res = (2*h/w)*res;
end

%b-spline pulse function
%INPUTS:
%t: current time
%w: width of pulse (starts at t=0, ends at t=h)
%h: height of pulse
%OUTPUTS:
%res: pulse evaluated at t
function res = b_spline_pulse(t,w,h)
    t = 4*t/w;
    b3 = (0<=t).*(t<1).*(t.^3)/4;
    t = t-1;
    b2 = (0<=t).*(t<1).*(-3*t.^3+3*t.^2+3*t+1)/4;
    t = t-1;
    b1 = (0<=t).*(t<1).*(3*t.^3-6*t.^2+4)/4;
    t = t-1;
    b0 = (0<=t).*(t<1).*(-t.^3+3*t.^2-3*t+1)/4;
    res = h*(b0+b1+b2+b3);
end
%b-spline pulse function (derivative)
%INPUTS:
%t: current time
%w: width of pulse (starts at t=0, ends at t=h)
%h: height of pulse
%OUTPUTS:
%res: derivative of pulse evaluated at t
function res = b_spline_pulse_derivative(t,w,h)
    t = 4*t/w;
    b3 = (0<=t).*(t<1).*(3*t.^2)/4;
    t = t-1;
    b2 = (0<=t).*(t<1).*(-9*t.^2+6*t+3)/4;
    t = t-1;
    b1 = (0<=t).*(t<1).*(9*t.^2-12*t)/4;
    t = t-1;
    b0 = (0<=t).*(t<1).*(-3*t.^2+6*t-3)/4;
    res = (4*h/w)*(b0+b1+b2+b3);
end




