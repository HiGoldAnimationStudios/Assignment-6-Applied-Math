%DISCRETE

function string_simulation_template01()

    %mypath1 = 'C:\Users\jvidaurrazaga\Downloads\';
    mypath1 = 'C:\Users\llin\Downloads\';
    fname='string_1.mp4';
    input_fname = [mypath1,fname];
    writerObj = VideoWriter(input_fname);
    open(writerObj);

    fig1=figure(1);

    num_masses = 3; %your code here
    total_mass = 10; %your code here
    tension_force = 5; %your code here
    string_length = 20; %your code here
    damping_coeff = 0.0001; %your code here

    dx = string_length/(num_masses+1);

    string_params = struct();
    string_params.n = num_masses;
    string_params.M = total_mass;
    string_params.Tf = tension_force;
    string_params.L = string_length;
    string_params.c = damping_coeff;
    string_params.dx = dx;

    mode_index = 3;

    [M_mat,K_mat] = construct_2nd_order_matrices(string_params);
    [Ur_mat,lambda_mat] = eig(K_mat,M_mat); %take UR_mat and plot the mode shape 
    
    omega_n = sqrt(lambda_mat(mode_index,mode_index));

    amplitude_Uf = 1;%your code here
    omega_Uf = omega_n; %your code here

    %list of x points (including the two endpoints)
    %xlist = linspace(0,string_length,num_masses+2);

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

    % Fehlberg = struct();
    % Fehlberg.C = [0, 1/4, 3/8, 12/13, 1, 1/2];
    % Fehlberg.B = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55;...
    % 25/216, 0, 1408/2565, 2197/4104, -1/5, 0];
    % Fehlberg.A = [0,0,0,0,0,0;...
    % 1/4, 0,0,0,0,0;...
    % 3/32, 9/32, 0,0,0,0;...
    % 1932/2197, -7200/2197, 7296/2197, 0,0,0;...
    % 439/216, -8, 3680/513, -845/4104, 0,0;...
    % -8/27, 2, -3544/2565, 1859/4104, -11/40, 0];
    % h_ref = 0.1;
    % error_desired = 10^-12;
    % p=5;

    %load string_params into rate function
    my_rate_func = @(t_in,V_in) string_rate_func01(t_in,V_in,string_params);
    
    %initial conditions
    U0 = zeros(num_masses,1);     
    dUdt0 = zeros(num_masses,1);
    V0 = [U0;dUdt0];

    %tspan = [0,100];
   
    %run the integration
    %[t_list,V_list,~, ~, ~, ~] = explicit_RK_variable_step_integration(my_rate_func,tspan,V0,h_ref,Fehlberg, p, error_desired);
    tspan = linspace(0,(20 + 0.75)*(2*pi)/omega_Uf,5000+1);
    [t_list,V_list] = ode45(my_rate_func,tspan,V0);
    V_list = V_list';

    hold on;

    maxU = max(max(V_list(1:num_masses,:)));
    minU = min(min(V_list(1:num_masses,:)));

    h = max([abs(maxU), abs(minU), amplitude_Uf]);

    %generating the mode shape plot
    x_distances = linspace(0, string_length, num_masses + 2);
    n_mode = Ur_mat(:,mode_index);
    mode_shape = [0; n_mode; 0];

    scale_factor = max(abs(maxU), abs(minU))/max(abs(mode_shape));
    mode_shape = mode_shape * scale_factor;

    plot(x_distances, mode_shape, 'o-', 'MarkerFaceColor', 'auto')

    ball_plot_struct = initialize_balls_plot();

    axis([0,string_length,-1.1*h,1.1*h]);
    xlabel('x-distance')
    ylabel('y-distance')
    legend('mode prediction', 'mode simulation')
    for k = 1:length(t_list)
        update_balls_plot(ball_plot_struct,V_list(:,k),t_list(k),string_params);
        drawnow;
        current_frame = getframe(fig1);
        %write the frame to the video
        writeVideo(writerObj,current_frame);
    end
   
    close(writerObj)
    
end