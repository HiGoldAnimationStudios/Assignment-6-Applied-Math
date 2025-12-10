function orion()
    num_masses = 5;
    total_mass = 10;
    tension_force = 2;
    string_length = 7;
    damping_coeff = 0.0001;
    dx = string_length/(num_masses+1);
    amplitude_Uf = 3;
    omega_Uf = 5;
    %list of x po5ints (including the two endpoints)
    xlist = linspace(0,string_length,num_masses+2);
    
    %generate the struct
    string_params = struct();


    string_params.n = num_masses;
    string_params.M = total_mass;
    
    string_params.Tf = tension_force;
    string_params.L = string_length;
    string_params.c = damping_coeff;
    string_params.dx = dx;

    mode_index=2;

    Uf_func = @(t_in) 0;
    dUfdt_func = @(t_in) 0;

    string_params.Uf_func = Uf_func;
    string_params.dUfdt_func = dUfdt_func;

    %load string_params into rate function
    my_rate_func = @(t_in,V_in) string_rate_func01(t_in,V_in,string_params);
    %initial conditions
    U0 = zeros(num_masses,1);
    dUdt0 = zeros(num_masses,1);
    V0 = [U0;dUdt0];
    %tspan = [0 5];
    %run the integration
    % [tlist,Vlist] = your_integrator(my_rate_func,tspan,V0,...);
    %your code to generate an animation of the system

    tlist=linspace(0,20*2*pi/omega_Uf,10000+1);

    Vresult=ode45(my_rate_func,tlist,V0);

    animate_string_with_pulse(tlist,Vresult,string_params)

    




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



function animate_string_with_pulse(tlist,Vresult,string_params)
    n=string_params.n;
    L=string_params.L;
    string_plot=plot(0,0,'o-','color','k','LineWidth',2, 'MarkerFaceColor','r',MarkerEdgeColor='r');
    
    maxU=max(max(Vresult(:,1:n)));
    
    minU=min(min(Vresult(:,1:n)));

    h=max(abs(maxU),abs(minU),string_params.Uf_amplitude);

    axis([0,L,-h,h])

    hold on
    xlabel("x")
    ylabel("U(t,x)")
    title('plot of vibrating string')


    for k=1:length(tlist)
        t=tlist(k);
        U=Vresult(k,1:n);
        Uf=string_params.Uf_func;
        U_padded=[0,U,Uf];
        x_list=linspace(0,L,n+2);
        drawnow;
    end
end