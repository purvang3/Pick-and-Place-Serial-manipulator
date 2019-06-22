#!/usr/bin/env python

# Copyright (C) 2017 Electric Movement Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *


def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:      
        ### Your FK code here
        q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8')
        d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
        a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
        alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')

        d01 = 0.75 # meters
        d34 = 1.5 # meters
        d4g = 0.303

        a12 = 0.35
        a23 = 1.25
        a34 = -0.054

        alpha12 = -pi/2
        alpha34 = -pi/2
        alpha45 = pi/2
        alpha56 = -pi/2

        # Create Modified DH parameters
        s = {alpha0:      0,   a0:   0, d1: d01, 
             alpha1: alpha12,  a1: a12, d2:   0,  q2:q2-pi/2,  
             alpha2:       0,  a2: a23, d3:   0,
             alpha3: alpha34,  a3: a34, d4: d34,
             alpha4: alpha45,  a4:   0, d5:   0,
             alpha5: alpha56,  a5:   0, d6:   0,
             alpha6:       0,  a6:   0, d7: d4g,   q7: 0}

        # Define Modified DH Transformation matrix
        def DH_Matrix (alpha,a,d,q):
            T_M = Matrix([[           cos(q),            -sin(q),           0,             a],
                          [ sin(q)*cos(alpha), cos(q)*cos(alpha), -sin(alpha), -sin(alpha)*d],
                          [ sin(q)*sin(alpha), cos(q)*sin(alpha),  cos(alpha),  cos(alpha)*d],
                          [                 0,                  0,           0,             1]])
            return T_M

        # Create individual transformation matrices
        T0_1 = DH_Matrix(alpha0,a0,d1,q1).subs(s)
        T1_2 = DH_Matrix(alpha1,a1,d2,q2).subs(s)
        T2_3 = DH_Matrix(alpha2,a2,d3,q3).subs(s)
        T3_4 = DH_Matrix(alpha3,a3,d4,q4).subs(s)
        T4_5 = DH_Matrix(alpha4,a4,d5,q5).subs(s)
        T5_6 = DH_Matrix(alpha5,a5,d6,q6).subs(s)
        T6_g = DH_Matrix(alpha6,a6,d7,q7).subs(s)

        T0_g = T0_1 * T1_2 * T2_3 * T3_4 * T4_5 * T5_6 * T6_g

        T0_2 = T0_1 * T1_2
        T0_3 = T0_2 * T2_3
        T0_4 = T0_3 * T3_4
        T0_5 = T0_4 * T4_5
        T0_6 = T0_5 * T5_6

        # Extract rotation matrices from the transformation matrices
        #m.extract([0,1,3],[0,1])
        #N= M[0:2,0:2]
        R0_1 = T0_1[0:3,0:3]
        R0_2 = T0_2[0:3,0:3]
        R0_3 = T0_3[0:3,0:3]
        R0_4 = T0_4[0:3,0:3]
        R0_5 = T0_5[0:3,0:3]
        R0_6 = T0_6[0:3,0:3]
        R0_g = T0_g[0:3,0:3]
        
        # Initialize service response
        joint_trajectory_list = []
        for x in xrange(0, len(req.poses)):
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()

            # Extract end-effector position and orientation from request
            # px,py,pz = end-effector position
            # roll, pitch, yaw = end-effector orientation
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z

            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, req.poses[x].orientation.y,
                req.poses[x].orientation.z, req.poses[x].orientation.w])
         
            ### Your IK code here 
            Rrpy =  Matrix([[cos(yaw)*cos(pitch), cos(yaw)*sin(pitch)*sin(roll)-sin(yaw)*cos(roll), cos(yaw)*sin(pitch)*cos(roll)+sin(yaw)*sin(roll)],
                            [sin(yaw)*cos(pitch), sin(yaw)*sin(pitch)*sin(roll)+cos(yaw)*cos(roll), sin(yaw)*sin(pitch)*cos(roll)-cos(yaw)*sin(roll)],
                            [-sin(pitch),         cos(pitch)*sin(roll),                             cos(pitch)*cos(roll)]])

            R_corr = Matrix([[0,0,1],
                            [0,-1,0],
                            [1,0,0]])

            Rrpy0_6 = Rrpy*(R_corr.T)

            EE=Matrix([[px],
                       [py],
                       [pz]])

            WC = EE - d4g * Rrpy0_6[:,2]
            

            theta1 = atan2(WC[1],WC[0])

            len1_6 = sqrt(pow((sqrt(WC[0] * WC[0]+WC[1] * WC[1]) - 0.35),2)+pow((WC[2]-0.75),2))
            len1_3 = a23
            len3_6 = sqrt(d34**2+a34**2)

            angle_temp1 = acos((len1_6 * len1_6 + len1_3 * len1_3 -len3_6 * len3_6) / (2 * len1_6 * len1_3))
            angle_temp2 = acos((len3_6 * len3_6 + len1_3 * len1_3 -len1_6 * len1_6) / (2 * len3_6 * len1_3))
            angle_temp3 = acos((len1_6 * len1_6 + len3_6 * len3_6 -len1_3 * len1_3) / (2 * len1_6 * len3_6))

            theta2 = pi/2 - angle_temp1 - atan2((WC[2]-0.75),(sqrt(WC[0] * WC[0]+WC[1] * WC[1]) - 0.35))
            theta3 = pi/2 - angle_temp2 - 0.036

            #R0_3 = T0_1[0:3,0:3] * T1_2[0:3,0:3] * T2_3[0:3,0:3]
            Rrpy0_3 = R0_3.evalf(subs={q1: theta1, q2: theta2, q3: theta3})
            #Rrpy3_6 = Rrpy0_3.inv("LU") * Rrpy0_6
            Rrpy3_6 = Rrpy0_3.T * Rrpy0_6

            #T3_6 = T3_4*T4_5*T5_6
            #R3_6 = T3_6[0:3,0:3]
            #print Rrpy3_6
            theta4 = atan2(Rrpy3_6[2,2], -Rrpy3_6[0,2])
            theta5 = atan2(sqrt(Rrpy3_6[0,2]*Rrpy3_6[0,2] + Rrpy3_6[2,2]*Rrpy3_6[2,2]),Rrpy3_6[1,2])
            theta6 = atan2(-Rrpy3_6[1,1],Rrpy3_6[1,0])
                
            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
            joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
            joint_trajectory_list.append(joint_trajectory_point)

        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()
