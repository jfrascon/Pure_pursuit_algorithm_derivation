clc;
clear all;
close all;

deltaT = 1/1000; % 1 ms to s

% Simulation time from t = 0 to t = tend
t = 0; 
tend = 60; % 1 min is 60 secs.

distance_to_wheel = 0.5; % 50 cm
rb = 1.4*distance_to_wheel; % robot's base radius

% I want to explore how to the robot behaves for diferent pairs of 
% (v_m_s, w_rad_s)
v_m_s = 10; %m/s
R = 0.7*distance_to_wheel;
w_rad_s = v_m_s/R;

% Initial pose for the robot.
robot_pose2D = [0; 0; 60/180*pi];

plot_robot(robot_pose2D, distance_to_wheel, rb, v_m_s, w_rad_s);
grid on;
axis([-5 5 -1.5 1.5]);
% Aspect ratio
daspect([1 1 1])
pbaspect([1 1 1])

% Container to store all positions covered by the robot.
poses2D = [robot_pose2D];

while t < tend    
    robot_pose2D = diff_kinematics(robot_pose2D, v_m_s, w_rad_s, deltaT); % current position
    poses2D = [poses2D robot_pose2D]; % Store position
    
    plot_robot(robot_pose2D, distance_to_wheel, rb, v_m_s, w_rad_s);
    hold on;
    grid on;
    plot(poses2D(1,:), poses2D(2,:), 'k--')
    axis([-5 5 -1.5 1.5]);
    % Aspect ratio
    daspect([1 1 1]);
    pbaspect([1 1 1]);
    pause(0.05);
    
    clf; % clear all the images in the current figure.
    
    t = t + deltaT; % update simulation time
end