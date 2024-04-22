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

#define pi 3.141592

const double t_hold = 0.5; 
const double t_swing1 = 1;
const double t_swing2 = 1;



char datapath[] = "data.csv";
char filename[] = "double_pendulum.xml";
char datafile[] = "data.csv";

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
const int data_frequency = 10; // frequency at which data is written to a file

// mouse interaction
bool button_left = false;
bool button_middle = false;
bool button_right = false;
double lastx = 0;
double lasty = 0;



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
    fprintf(fid, "PE, KE, TE, ");
    fprintf(fid, "q1, q2, ");

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
    fprintf(fid, "%f, %f, %f, ", d->energy[0], d->energy[1], d->energy[0] + d->energy[1]);
    fprintf(fid, "%f, %f ", d->qpos[0], d->qpos[1]);
    //Don't remove the newline
    fprintf(fid, "\n");
}

//void init_controller(const mjModel* m, mjData* D)
//{
 //fsm_state = fsm_hold;   
//}
double tustin_derivative(double input, double input_old, double output_old, double cutoff_freq)
{
    double time_const = 1 / (2 * pi * cutoff_freq);
    double output = 0;

    output = (2 * (input - input_old) - (0.0001 - 2 * time_const) * output_old) / (0.0001 + 2 * time_const);

    return output;
}
// 모터 위치제어에 필요한 2차 low pass filter * inertia
double tustin_2nd_derivative(double input_old2, double output_old1, double output_old2, double J,double cutoff_freq, double damping)
{
    double wc= 2*pi*cutoff_freq;
    double Q = 1/(2*damping);
    double output = 0;
    output = J*input_old2-1/(Q*wc)*output_old1-1/(pow(wc,2))*output_old2;
    return output;
}
