# RobotKinematics
This MATLAB script is a function that calculates the inverse kinematics of a standing series robotic arm with 5 joints (UR5) using Newton-Raphson algorithm and prints out the iterations. For the UR5 robot, The home configuration along with the desired configuration of the end-effector M in form of transformation matrices, the velocity vectors (angular(omega) and linear (V)) in form of a screw axes in the end-effector frame, the absolute error and the initial guess must all be given as user input. The most notable output of the function is the joint values row vector sorted in a matrix with iterations as rows. The function can be customized as desired. If you need further assistance, contact me. 

The animation video records the simulation of the function by the V-REP, a robot simulation platform. 


Made by ahmed11406
https://github.com/ahmed11406/RobotKinematics
