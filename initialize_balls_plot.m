function ball_plot_struct = initialize_balls_plot(string_params)
    
    num_count = string_params.n;
    xpoints = 0 : string_params.dx : string_params.dx*(num_count - 1);

    for i = 1:length(xpoints)
        ball_plot_struct = struct();
        ball_plot_struct.line_plot = plot(xpoints(i),0,'k','linewidth',2);
        ball_plot_struct.point_plot = plot(xpoints(i),0,'ro','markerfacecolor','r','markersize',7);
    
    end
end