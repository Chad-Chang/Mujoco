#pragma once
#include<stdbool.h> //for bool
//#include<unistd.h> //for usleep
//#include <math.h>

#include <GLFW/glfw3.h>
#include <mujoco/mujoco.h>
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "PID.h" 
#define ndof 2
#define fsm_hold 0
#define fsm_swing1 1
#define fsm_swing2 2
#define fsm_stop 3
#define G 9.81

#define pi 3.141592

const double t_hold = 0.5; 
const double t_swing1 = 1;
const double t_swing2 = 1;
double err_old=  0; double err = 0;  
double ref = -pi/2;
double u_c = 0;
double dense_M[ndof * ndof] = { 0 };
double d_hat[2] = {0}; // disturbance
double theta[3] = {0}; // output angle from encoder
double theta_acc[3] = {0};
double u_d[3] = {0}; // uc_- u_d
double theta_vel[3] = {0};
double dist_freq = 1;
double LPF_freq = 30;
double perturb; // disturbance
float amplitude_perturb = 5;
bool dob_switch = 1; 
double ext_force[2] = {0}; // external force(x,y)
double ext_force_hat[2] = {0};
double l_force = 1; // 힘이 작용하는 길이 (x^2 + y^2 = l^2) 
double est_torque[3] ={0};
double admi = 0;

char filename[] = "double_pendulum.xml";
char datafile[] = "data/DOB.csv";

// MuJoCo data structures
mjModel* m = NULL;                  // MuJoCo model
mjData* d = NULL;                   // MuJoCo data
mjvCamera cam;                      // abstract camera
mjvOption opt;                      // visualization options
mjvScene scn;                       // abstract scene
mjrContext con;                     // custom GPU context

// Simulation End Time
double simend = 5;
// const int nv = 2; // DoF of system 

// Data Writing
FILE* fid;
int loop_index = 0;
const int data_frequency = 4; // frequency at which data is written to a file

// mouse interaction
bool button_left = false;
bool button_middle = false;
bool button_right = false;
double lastx = 0;
double lasty = 0;
double Ts = 0.0001;
double cutoff = 150;

// keyboard callback
void keyboard(GLFWwindow* window, int key, int scancode, int act, int mods)
{
    // backspace: reset simulation
    if (act == GLFW_PRESS && key == GLFW_KEY_BACKSPACE)
    {
        mj_resetData(m, d);
        mj_forward(m, d);
    }
}
// mouse button callback
void mouse_button(GLFWwindow* window, int button, int act, int mods)
{
    // update button state
    button_left = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS);
    button_middle = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS);
    button_right = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS);

    // update mouse position
    glfwGetCursorPos(window, &lastx, &lasty);
}
// mouse move callback
void mouse_move(GLFWwindow* window, double xpos, double ypos)
{
    // no buttons down: nothing to do
    if (!button_left && !button_middle && !button_right)
        return;

    // compute mouse displacement, save
    double dx = xpos - lastx;
    double dy = ypos - lasty;
    lastx = xpos;
    lasty = ypos;

    // get current window size
    int width, height;
    glfwGetWindowSize(window, &width, &height);

    // get shift key state
    bool mod_shift = (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
        glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS);

    // determine action based on mouse button
    mjtMouse action;
    if (button_right)
        action = mod_shift ? mjMOUSE_MOVE_H : mjMOUSE_MOVE_V;
    else if (button_left)
        action = mod_shift ? mjMOUSE_ROTATE_H : mjMOUSE_ROTATE_V;
    else
        action = mjMOUSE_ZOOM;

    // move camera
    mjv_moveCamera(m, action, dx / height, dy / height, &scn, &cam);
}
// scroll callback
void scroll(GLFWwindow* window, double xoffset, double yoffset)
{
    // emulate vertical mouse motion = 5% of window height
    mjv_moveCamera(m, mjMOUSE_ZOOM, 0, -0.05 * yoffset, &scn, &cam);
}


// holders of one step history of time and position to calculate dertivatives
mjtNum position_history = 0;
mjtNum previous_time = 0;

// controller related variables
float_t ctrl_update_freq = 100;
mjtNum last_update = 0.0;
mjtNum ctrl;

// Position Servo
void set_position_servo(const mjModel* m, int actuator_no, double kp) {
    m->actuator_gainprm[10 * actuator_no + 0] = kp;
    m->actuator_biasprm[10 * actuator_no + 1] = -kp;
}

// Velocity Servo
void set_velocity_servo(const mjModel* m, int actuator_no, double kv) {
    m->actuator_gainprm[10 * actuator_no + 0] = kv;
    m->actuator_biasprm[10 * actuator_no + 2] = -kv;
}

// Torque Control
void set_torque_control(const mjModel* m, int actuator_no, int flag) {
    if (flag == 0)    // off
        m->actuator_gainprm[10 * actuator_no + 0] = 0;
    else
        m->actuator_gainprm[10 * actuator_no + 0] = 1;
}

//****************************
//This function is called once and is used to get the headers
void init_save_data() // csv파일의 데이터 명을 지정하는 함수 -> 한번만 실행
{
    //write name of the variable here (header)
    fprintf(fid, "t, ");
    fprintf(fid, "dist, e_dist, dist_err, "); // disturbace error
    fprintf(fid, "input, output, error, gravity compensation,"); // reference error
    fprintf(fid,  "ext_force_x_hat, ext_force_y_hat, torque_hat, ref, admittance,");

    //Don't remove the newline
    fprintf(fid, "\n");
}

////***************************
////This function is called at a set frequency, put data here
void save_data(const mjModel* m, mjData* d) 
{   
    //data here should correspond to headers in init_save_data()
    //seperate data by a space %f followed by space
    fprintf(fid, "%f, ", d->time);
    fprintf(fid, "%f, %f, %f, ", perturb, d_hat[0], perturb- d_hat[0]);
    fprintf(fid, "%f, %f, %f, %f, ", ref,d-> qpos[0], ref - d-> qpos[0],0.5*G*sin(d->qpos[0]));
    fprintf(fid, "%f,%f,%f,%f, %f", ext_force_hat[0], ext_force_hat[1], est_torque[0], ref, admi);

    //Don't remove the newline
    fprintf(fid, "\n");
}

double tustin_derivative(double input, double input_old, double output_old, double cutoff_freq)
{
    double time_const = 1 / (2 * pi * cutoff_freq);
    double output = 0;

    output = (2 * (input - input_old) - (0.0001 - 2 * time_const) * output_old) / (0.0001 + 2 * time_const);

    return output;
}

double t_k =0 ;
double t_kold= 0 ;
int contol_loop = 4; // 0.0001ms초를 보장해주는 놈
void loop_tcheck()// loop time check해주는 놈
{ 
    t_k = d->time;
    printf("%f \n", t_k-t_kold);
    t_kold = t_k;
}

// 1 / (tau*s+1)
double lowpassfilter(double input, double input_old, double output_old, double cutoff_freq) 
{
    double time_const = 1 / (2 * pi * cutoff_freq);
    double output = 0;

    output = (Ts * (input + input_old) - (Ts - 2 * time_const) * output_old) / (Ts + 2 * time_const);

    return output;
}