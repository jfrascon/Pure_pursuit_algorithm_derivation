function [pose2D] = diff_kinematics(pose2D, v_m_s, w_rad_s, deltaT)

  pose2D = pose2D + [v_m_s*cos(pose2D(3))*deltaT; v_m_s*sin(pose2D(3))*deltaT; w_rad_s*deltaT];

end
