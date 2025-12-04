function string_simulation_template01()
    num_masses = 2; %your code here
    total_mass = 10; %your code here
    tension_force = 5; %your code here
    string_length = 7; %your code here
    damping_coeff = 0.67; %your code here
    
    dx = string_length/(num_masses+1);

    amplitude_Uf = 1;%your code here
    omega_Uf = 1; %your code here

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

    Fehlberg = struct();
    Fehlberg.C = [0, 1/4, 3/8, 12/13, 1, 1/2];
    Fehlberg.B = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55;...
    25/216, 0, 1408/2565, 2197/4104, -1/5, 0];
    Fehlberg.A = [0,0,0,0,0,0;...
    1/4, 0,0,0,0,0;...
    3/32, 9/32, 0,0,0,0;...
    1932/2197, -7200/2197, 7296/2197, 0,0,0;...
    439/216, -8, 3680/513, -845/4104, 0,0;...
    -8/27, 2, -3544/2565, 1859/4104, -11/40, 0];
    h_ref = 2.8;
    error_desired = 10^-6;
    p=5;

    %load string_params into rate function
    my_rate_func = @(t_in,V_in) string_rate_func01(t_in,V_in,string_params);
    
    %initial conditions
    U0 = zeros(num_masses,1);     %your code here
    dUdt0 = zeros(num_masses,1);%your code here
    V0 = [U0;dUdt0];

    tspan = [0,5];%your code here
   
    %run the integration
    [t_list,V_list,~, ~, ~, ~] = explicit_RK_variable_step_integration(my_rate_func,tspan,V0,h_ref,Fehlberg, p, error_desired);
    
    ball_plot_struct = initialize_balls_plot;

    N = length(t_list);

    for k = 1:N
        update_balls_plot(ball_plot_struct,V_list(k))
    end
end