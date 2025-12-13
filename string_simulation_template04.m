%harmonics?

function string_simulation_template04() 
    num_masses = 50; %your code here
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
    
    mode_index=3; 

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
    
    U0 = zeros(num_masses,1);     
    dUdt0 = zeros(num_masses,1);
    V0 = [U0;dUdt0];

    my_rate_func = @(t_in,V_in) string_rate_func01(t_in,V_in,string_params);

    % tspan = linspace(0,3*string_length/c,501);
    tspan = linspace(0,20*(2*pi/omega_Uf),2001);
    [t_list,V_list] = ode45(my_rate_func,tspan,V0);
    V_list = V_list';
    
    q = max(max(abs(V_list(1:num_masses,:))));
    
    n=5;
    n=50;
    mode_list=linspace(1,n,n);
    omega_list_discrete=linspace(1,n,n);
    omega_list_continuous=linspace(1,n,n);
    

    figure(1);
    for j=1:length(mode_list)
        mode_index=mode_list(j);
        string_params.n=n;
        string_params.dx=string_length/(n+1);

        [M_mat,K_mat] = construct_2nd_order_matrices(string_params);
        [Ur_mat,lambda_mat] = eig(K_mat,M_mat); 

        omega_list_discrete(j)=sqrt(lambda_mat(mode_index,mode_index));
        omega_n_spatial=pi*mode_index/string_length;
        omega_list_continuous(j)=c*omega_n_spatial;
    end

    hold on
    title("Mode Index vs Frequency (50 masses)")
    xlabel("mode index")
    ylabel("omega")
    legend();
    plot(mode_list,omega_list_discrete,"ro-","DisplayName","Discrete")
    plot(mode_list,omega_list_continuous,"bo-","DisplayName","Continuous")
end

