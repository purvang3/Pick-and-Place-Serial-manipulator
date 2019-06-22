from sympy import *
from time import time
from mpmath import radians
import tf

'''
Format of test case is [ [[EE position],[EE orientation as quaternions]],[WC location],[joint angles]]
You can generate additional test cases by setting up your kuka project and running `$ roslaunch kuka_arm forward_kinematics.launch`
From here you can adjust the joint angles to find thetas, use the gripper to extract positions and orientation (in quaternion xyzw) and lastly use link 5
to find the position of the wrist center. These newly generated test cases can be added to the test_cases dictionary.
'''

test_cases = {1:[[[2.16135,-1.42635,1.55109],
                  [0.708611,0.186356,-0.157931,0.661967]],
                  [1.89451,-1.44302,1.69366],
                  [-0.65,0.45,-0.36,0.95,0.79,0.49]],
              2:[[[-0.56754,0.93663,3.0038],
                  [0.62073, 0.48318,0.38759,0.480629]],
                  [-0.638,0.64198,2.9988],
                  [-0.79,-0.11,-2.33,1.94,1.14,-3.68]],
              3:[[[-1.3863,0.02074,0.90986],
                  [0.01735,-0.2179,0.9025,0.371016]],
                  [-1.1669,-0.17989,0.85137],
                  [-2.99,-0.12,0.94,4.06,1.29,-4.12]],
              4:[],
              5:[]}


def test_code(test_case):
    ## Set up code
    ## Do not modify!
    x = 0
    class Position:
        def __init__(self,EE_pos):
            self.x = EE_pos[0]
            self.y = EE_pos[1]
            self.z = EE_pos[2]
    class Orientation:
        def __init__(self,EE_ori):
            self.x = EE_ori[0]
            self.y = EE_ori[1]
            self.z = EE_ori[2]
            self.w = EE_ori[3]

    position = Position(test_case[0][0])
    orientation = Orientation(test_case[0][1])

    class Combine:
        def __init__(self,position,orientation):
            self.position = position
            self.orientation = orientation

    comb = Combine(position,orientation)

    class Pose:
        def __init__(self,comb):
            self.poses = [comb]

    req = Pose(comb)
    start_time = time()
    
    ########################################################################################
    ## 

    ## Insert IK code here!
    
    ### Your FK code here
    # Create symbols
    q1, q2, q3, q4, q5, q6, qg = symbols('q1:8')
    d1, d2, d3, d4, d5, d6, dg = symbols('d1:8')
    a0, a1, a2, a3, a4, a5, ag = symbols('a0:7')
    alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alphag = symbols('alpha0:7')

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
         alphag:       0,  ag:   0, dg: d4g,   qg: 0}

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
    T6_g = DH_Matrix(alphag,ag,dg,qg).subs(s)

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
    print Rrpy3_6
    theta4 = atan2(Rrpy3_6[2,2], -Rrpy3_6[0,2])
    theta5 = atan2(sqrt(Rrpy3_6[0,2]*Rrpy3_6[0,2] + Rrpy3_6[2,2]*Rrpy3_6[2,2]),Rrpy3_6[1,2])
    theta6 = atan2(-Rrpy3_6[1,1],Rrpy3_6[1,0])

    print ("\ntheta1 = %f" %theta1)
    print ("\ntheta2 = %f" %theta2)
    print ("\ntheta3 = %f" %theta3)
    print ("\ntheta4 = %f" %theta4)
    print ("\ntheta5 = %f" %theta5)
    print ("\ntheta6 = %f" %theta6)

    #theta4 = atan2(Rrpy3_6[2,2], -Rrpy3_6[0,2])
    #theta5 = atan2(sqrt(Rrpy3_6[0,2]*Rrpy3_6[0,2]+Rrpy3_6[2,2]*Rrpy3_6[2,2]),Rrpy3_6[1,2])
    #theta6 = atan2(-Rrpy3_6[1,1],Rrpy3_6[1,0])
    ## 
    ########################################################################################
    
    ########################################################################################
    ## For additional debugging add your forward kinematics here. Use your previously calculated thetas
    ## as the input and output the position of your end effector as your_ee = [x,y,z]

    ## (OPTIONAL) YOUR CODE HERE!
    #-0.65,0.45,-0.36,0.95,0.79,0.49
    # theta1 = -0.65
    # theta2 = 0.45
    # theta3 = -0.36
    # theta4 = 0.95
    # theta5 = 0.79
    # theta6 = 0.49

    FK = T0_g.evalf(subs={q1: theta1, q2: theta2, q3: theta3, q4: theta4, q5: theta5, q6: theta6})

    ## End your code input for forward kinematics here!
    ########################################################################################

    ## For error analysis please set the following variables of your WC location and EE location in the format of [x,y,z]
    your_wc = [WC[0],WC[1],WC[2]] # <--- Load your calculated WC values in this array
    your_ee = [FK[0,3], FK[1,3], FK[2,3]] # <--- Load your calculated end effector value from your forward kinematics
    ########################################################################################

    ## Error analysis
    print ("\nTotal run time to calculate joint angles from pose is %04.4f seconds" % (time()-start_time))

    # Find WC error
    if not(sum(your_wc)==3):
        wc_x_e = abs(your_wc[0]-test_case[1][0])
        wc_y_e = abs(your_wc[1]-test_case[1][1])
        wc_z_e = abs(your_wc[2]-test_case[1][2])
        wc_offset = sqrt(wc_x_e**2 + wc_y_e**2 + wc_z_e**2)
        print ("\nWrist error for x position is: %04.8f" % wc_x_e)
        print ("Wrist error for y position is: %04.8f" % wc_y_e)
        print ("Wrist error for z position is: %04.8f" % wc_z_e)
        print ("Overall wrist offset is: %04.8f units" % wc_offset)

    # Find theta errors
    t_1_e = abs(theta1-test_case[2][0])
    t_2_e = abs(theta2-test_case[2][1])
    t_3_e = abs(theta3-test_case[2][2])
    t_4_e = abs(theta4-test_case[2][3])
    t_5_e = abs(theta5-test_case[2][4])
    t_6_e = abs(theta6-test_case[2][5])
    print ("\nTheta 1 error is: %04.8f" % t_1_e)
    print ("Theta 2 error is: %04.8f" % t_2_e)
    print ("Theta 3 error is: %04.8f" % t_3_e)
    print ("Theta 4 error is: %04.8f" % t_4_e)
    print ("Theta 5 error is: %04.8f" % t_5_e)
    print ("Theta 6 error is: %04.8f" % t_6_e)
    print ("\n**These theta errors may not be a correct representation of your code, due to the fact \
           \nthat the arm can have muliple positions. It is best to add your forward kinmeatics to \
           \nconfirm whether your code is working or not**")
    print (" ")

    # Find FK EE error
    if not(sum(your_ee)==3):
        ee_x_e = abs(your_ee[0]-test_case[0][0][0])
        ee_y_e = abs(your_ee[1]-test_case[0][0][1])
        ee_z_e = abs(your_ee[2]-test_case[0][0][2])
        ee_offset = sqrt(ee_x_e**2 + ee_y_e**2 + ee_z_e**2)
        print ("\nEnd effector error for x position is: %04.8f" % ee_x_e)
        print ("End effector error for y position is: %04.8f" % ee_y_e)
        print ("End effector error for z position is: %04.8f" % ee_z_e)
        print ("Overall end effector offset is: %04.8f units \n" % ee_offset)




if __name__ == "__main__":
    # Change test case number for different scenarios
    test_case_number = 3

    test_code(test_cases[test_case_number])
