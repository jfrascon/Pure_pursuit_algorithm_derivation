function plot_robot(pose2D, distance_to_wheel, rb, v_m_s, w_rad_s)

% Plots a circular robot with its body-frame. body frame's x-axis is red,
% y-axis is green.

% pose2D: robot's pose in a 2D world
% distance_to_wheel: distance from the robot's center to each wheel
% rb: robot's base radius.

length_arrow = 0.5*distance_to_wheel;

angles = 0:0.05:2*pi;
footprint = pose2D(1:2)*ones(size(angles)) + rb*[cos(angles); sin(angles)];

% Rotation matrix around Z-axis
R = [cos(pose2D(3)) -sin(pose2D(3)); sin(pose2D(3)) cos(pose2D(3))];

% In local coordinates.
posfin_xarrow_local = [length_arrow; 0];
posfin_yarrow_local = [0; length_arrow;];

% Cross repressenting the left wheel, in local coordinates, so 4 points
% (two lines intersecting at 90 degrees)
posini_lwheel_vert_local = [0; distance_to_wheel - 0.1*distance_to_wheel];
posend_lwheel_vert_local = [0; distance_to_wheel + 0.1*distance_to_wheel];
posini_lwheel_hor_local = [-0.1*distance_to_wheel; distance_to_wheel];
posend_lwheel_hor_local = [+0.1*distance_to_wheel; distance_to_wheel];

pos_lwheel_local = [posini_lwheel_vert_local posend_lwheel_vert_local posini_lwheel_hor_local posend_lwheel_hor_local];

% Cross repressenting the left wheel, in local coordinates, so 4 points
% (two lines intersecting at 90 degrees)
posini_rwheel_vert_local = [0; -distance_to_wheel + 0.1*distance_to_wheel];
posend_rwheel_vert_local = [0; -distance_to_wheel - 0.1*distance_to_wheel];
posini_rwheel_hor_local = [-0.1*distance_to_wheel; -distance_to_wheel];
posend_rwheel_hor_local = [+0.1*distance_to_wheel; -distance_to_wheel];

pos_rwheel_local = [posini_rwheel_vert_local posend_rwheel_vert_local posini_rwheel_hor_local posend_rwheel_hor_local];

% Coordinates of the turning radius in local
turning_radius_local = [0; v_m_s/w_rad_s];

% Transform to global coordinates.
posfin_xarrow = R * posfin_xarrow_local + pose2D(1:2);
posfin_yarrow = R * posfin_yarrow_local + pose2D(1:2);
pos_lwheel = R * pos_lwheel_local + pose2D(1:2).*ones(2,4);
pos_rwheel = R * pos_rwheel_local + pose2D(1:2).*ones(2,4);
turning_radius = R * turning_radius_local + pose2D(1:2);

% Plots
plot(footprint(1,:), footprint(2,:))
hold on;
arrow(pose2D(1:2), posfin_xarrow, 30, 'BaseAngle', 60,'EdgeColor', 'r', 'FaceColor', 'r');
arrow(pose2D(1:2), posfin_yarrow, 30, 'BaseAngle', 60,'EdgeColor', 'g', 'FaceColor', 'g');
line([pos_lwheel(1,1) pos_lwheel(1,2)], [pos_lwheel(2,1) pos_lwheel(2,2)])
line([pos_lwheel(1,3) pos_lwheel(1,4)], [pos_lwheel(2,3) pos_lwheel(2,4)])
line([pos_rwheel(1,1) pos_rwheel(1,2)], [pos_rwheel(2,1) pos_rwheel(2,2)])
line([pos_rwheel(1,3) pos_rwheel(1,4)], [pos_rwheel(2,3) pos_rwheel(2,4)])
plot(turning_radius(1), turning_radius(2), '.', 'Markersize' ,36);
end