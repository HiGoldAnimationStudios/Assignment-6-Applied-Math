function ball_plot_struct = initialize_balls_plot()

    ball_plot_struct = struct();
    ball_plot_struct.line_plot = plot(0,0,'k','linewidth',2);
    ball_plot_struct.point_plot = plot(0,0,'ro','markerfacecolor','r','markersize',7);
    
end