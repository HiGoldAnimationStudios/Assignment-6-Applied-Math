function update_balls_plot(ball_plot_struct,V)
    l=length(V)/2;
    U = V(1:l);
    
    set(ball_plot_struct.line_plot,'xdata',plot_pts(1,:),'ydata',plot_pts(2,:));
    set(ball_plot_struct.point_plot,'xdata',plot_pts(1,:),'ydata',plot_pts(2,:));
end