function update_balls_plot(ball_plot_struct,V,t,string_params)
    % % X points
    % x_points=[0];
    % for i=1:string_params.n
    %     x_points(end+1)=string_params.dx*i;
    % end
    % x_points(end+1)=string_params.L;


    x_points = linspace(0,string_params.L,string_params.n+2);

    % Y points
    l=length(V)/2;
    U = V(1:l)';
    U = [0, U, string_params.Uf_func(t)];
   
    set(ball_plot_struct.line_plot,'xdata',x_points,'ydata',U);
    set(ball_plot_struct.point_plot,'xdata',x_points,'ydata',U);
end