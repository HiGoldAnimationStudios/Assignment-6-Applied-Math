function string_simulation_template01()
    num_masses = 2;%your code here
    total_mass = 10;%your code here
    tension_force = 5;%your code here
    string_length = 1;%your code here
    damping_coeff = 0.7;%your code here
    dx = string_length/(num_masses+1);
    amplitude_Uf = 1;%your code here
    omega_Uf = 1;%your code here
    %list of x points (including the two endpoints)
    xlist = linspace(0,string_length,num_masses+2);
    Uf_func = @(t_in) amplitude_Uf*cos(omega_Uf*t_in);
    dUfdt_func = @(t_in) -omega_Uf*amplitude_Uf*sin(omega_Uf*t_in);
    %generate the struct
    string_params = struct();
    string_params.n = num_masses;
    string_params.M = total_mass;
    string_params.Uf_func = Uf_func;
    string_params.dUfdt_func = dUfdt_func;
    string_params.Tf = tension_force;
    string_params.L = string_length;
    string_params.c = damping_coeff;
    string_params.dx = dx;
    %load string_params into rate function
    my_rate_func = @(t_in,V_in) string_rate_func01(t_in,V_in,string_params);
    %initial conditions
    U0 = 0;%your code here
    dUdt0 = 0;%your code here
    V0 = [U0;dUdt0];
    tspan = [0,5];%your code here
    %run the integration
    %[tlist,Vlist] = your_integrator(my_rate_func,tspan,V0,...);
    %your code to generate an animation of the system
    [t_list_nonlinear,V_list_nonlinear,~, ~, ~, ~] = explicit_RK_variable_step_integration(my_rate_func,tspan,V0,h_ref,Fehlberg, p, error_desired);
    %[t_list_linear,V_list_linear,~, ~, ~, ~] = explicit_RK_variable_step_integration(my_linear_rate,tspan,V0,h_ref,Fehlberg, p, error_desired);
    %[t_list,X_list] = ode45(my_rate_func,tspan,V0);

    delta_x_list_nonlinear=V_list_nonlinear(1,:)-Veq(1);
    delta_y_list_nonlinear=V_list_nonlinear(2,:)-Veq(2);
    delta_theta_list_nonlinear=V_list_nonlinear(3,:)-Veq(3);
    delta_x_list_linear=V_list_linear(1,:)-Veq(1);
    delta_y_list_linear=V_list_linear(2,:)-Veq(2);
    delta_theta_list_linear=V_list_linear(3,:)-Veq(3);

    
    figure(1)

    hold on
    plot(t_list_linear,delta_x_list_linear,"r", "DisplayName","\Deltax linear")
    plot(t_list_linear,delta_y_list_linear,"b", "DisplayName","\Deltay linear")
    plot(t_list_linear,delta_theta_list_linear,"g", "DisplayName","\Delta\theta linear")
    plot(t_list_nonlinear,delta_x_list_nonlinear,"r--", "DisplayName","\Deltax nonlinear")
    plot(t_list_nonlinear,delta_y_list_nonlinear,"b--", "DisplayName","\Deltay nonlinear")
    plot(t_list_nonlinear,delta_theta_list_nonlinear,"g--", "DisplayName","\Delta\theta nonlinear")
    legend("Location","best")
    xlim([0 5])
    xlabel("time")
    ylabel("displacement")
    title("linear vs non linear from tiny pertubation \epsilon=0.0067")
    

end