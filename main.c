// #pragma once

#include "funcs.h" // 사용자 정의 헤더파일 <>하면 안됨.
// namspace std
PID pid1;
PID pid2;
struct ParamModel_{
    double q[2];
    double qd[2];
    double qdd[2];
    double ctrl_input[2];
};
void state_init(const mjModel *m, mjData *d, ParamModel_ *model, int dof)
{
    for(int i =0; i<2 ; i++)
    {   
        model->q[i] = 0;
        model->qd[i] = 0;
        model->qdd[i] = 0;
        model->ctrl_input[i] = 0;
        //asdf
    }
}
void update_state(const mjModel *m, mjData *d, ParamModel_ *model, int dof)
{
    for(int i = 0 ; i < 1 ; i++ )
    {
        model->q[i] = d -> qpos[i];
        model->qd[i] = d -> qvel[i];
        model->qdd[i] = d -> qacc[i];
        // model->ctrl_input[i] = input[i];
    }
}

ParamModel_ dPendulum;
double alpha[3] = {0};
double beta[3] = {0};
bool fob = 1; 

double J_trans_inv[2] = {0}; 
double ctrl_torque[3] = {0};
double ctrl_torque_s1[3] =  {0};
double ctrl_torque_s2[3] =  {0};
double delta_pos[3] = {0}; 
double mass = 1; 
double nominal_plant = 0;

void admittance(double input, double input_old1, double input_old2, double output_old1, double output_old2, double K,double zeta, double wn)
{    
    // wn = sqrt(K/M), zeta = C / sqrt(4*M*K)
    double nominator = (pow(wn,2)/K) * (input_old2+2*input_old1 + input) - (4/pow(Ts,2)-4*wn*zeta/Ts + pow(wn,2))*(output_old2) - (2*pow(wn,2)-8/(pow(Ts,2)))*output_old1;  
    double denominator  = (4/pow(Ts,2) +4*wn*zeta/Ts + pow(Ts,2));
    admi = nominator / denominator;
    delta_pos[0] = admi;
    // printf("admittance = %f\n", admi);
}
void admittance2(double input, double input_old1, double input_old2, double output_old1, double output_old2, double J,double B, double K)
{   
    double zeta = B/sqrt(4*K*J);
    double wn = sqrt(K/J);
    double nominator = (1/J) * (input_old2+2*input_old1 + input) - (4/pow(Ts,2)-4*wn*zeta/Ts + pow(wn,2))*(output_old2) - (2*pow(wn,2)-8/(pow(Ts,2)))*output_old1;  
    double denominator  = 4/pow(Ts,2) +4*wn*zeta/Ts + pow(Ts,2);
    admi = nominator / denominator;
    delta_pos[0] = admi;
    printf("delta_pos = %f\n", admi);
}
void admittance3(double delta_force, double K)
{   
    admi = delta_force / K;
    delta_pos[0] = admi;
    printf("delta_pos = %f, %f\n", admi, delta_force);
}

void mycontroller(const mjModel* m, mjData* d)  // 제어주기 0.000025임
{   
    double u_g = 0.5*G*sin(d->qpos[0]);

        
    // delta_pos[0] = -pi/2 - d->qpos[0];
// single pendulum
    if(loop_index % contol_loop ==0) // sampling time 0.0001
    {   
        double M[ndof][ndof] = { 0 };   // 2x2 inertia matrix -> 알아서 계산해줌.
        mj_fullM(m, dense_M, d->qM);
        
        M[0][0] = dense_M[0];
        M[0][1] = dense_M[1];
        M[1][0] = dense_M[2];
        M[1][1] = dense_M[3];
        printf("%f,%f,%f,%f\n", M[0][0], M[0][1], M[1][0], M[1][1]);
        // nominal_plant = M[0][0];
        perturb = amplitude_perturb*sin(dist_freq*2*pi*d->time);
        // loop_tcheck();

        perturb = 0;

        ref = theta_des;

        err = ref - theta[0];
        pid1.set_gain(15,0,1);
            // control torque 
        u_c = pid1.compute_PID(err,err_old, Ts, cutoff);
        u_d[0] = u_c-d_hat[0];
        d->ctrl[0] = u_d[0] + perturb - u_g;// - u_g;//u_d[0] ;
        ctrl_torque[0]= u_d[0];

        ref = -pi/2;

        err_old = err;
        pid1.update_PID();

        theta[0] = d->qpos[0];
        theta_vel[0] = tustin_derivative(theta[0],theta[1], theta_vel[1], 300);
        theta_acc[0] = tustin_derivative(theta_vel[0],theta_vel[1], theta_acc[1], LPF_freq);
        alpha[0] = lowpassfilter(u_d[0], u_d[1], alpha[1], LPF_freq);

        

        if(dob_switch)
        {
            d_hat[0] = (0.25+0.5*M[0][0])*theta_acc[0]-beta[0];
            // printf("disturb = %f, gravity = %f, error = %f \n", d_hat[0], u_g, d_hat[0]- u_g);
        }

        if(fob)
        {
            J_trans_inv[0] = cos(d->qpos[0])/l_force;
            J_trans_inv[1] = -sin(d->qpos[0])/l_force;
            ctrl_torque_s1[0] = lowpassfilter(ctrl_torque[0],ctrl_torque[1],ctrl_torque_s1[1],LPF_freq);
            ext_force_hat[0] = J_trans_inv[0]*(ctrl_torque_s1[0]-(M[0][0])*theta_acc[0]);
            ext_force_hat[1] = J_trans_inv[1]*(ctrl_torque_s1[0]-(M[0][0])*theta_acc[0]);
            est_torque[0] = - ctrl_torque_s1[0] + M[0][0]*theta_acc[0];
            est_torque[0] = 30;// cos(2*pi**d->time);
            
            if(d->time >3) est_torque[0] = 0;
            admittance(est_torque[0], est_torque[1], est_torque[2],delta_pos[1], delta_pos[2], 100, 0.5, 100);  // 고유진동수, 
            
            // ref = ref + delta_pos[0];


            est_torque[2] = est_torque[1]; est_torque[1] = est_torque[0];
            delta_pos[2] = delta_pos[1]; delta_pos[1] = delta_pos[0]; 
        }


        theta[2] = theta[1]; theta[1] = theta[0];
        theta_vel[2] = theta_vel[1]; theta_vel[1] = theta_vel[0];
        theta_acc[2] = theta_acc[1]; theta_acc[1] = theta_acc[0];
        alpha[1] = alpha[0];
        beta[1] = beta[0];
        u_d[2] = u_d[1]; u_d[1] = u_d[0];
        ctrl_torque[1] = ctrl_torque[0]; 
        ctrl_torque_s1[1] = ctrl_torque_s1[0]; 
        ctrl_torque_s2[1] = ctrl_torque_s2[0]; 
            }

    if (loop_index % data_frequency == 0) {     // loop_index를 data_frequency로 나눈 나머지가 0이면 데이터를 저장.
        save_data(m, d);
        
    }
    //loop_index += 1;
    loop_index = loop_index + 1;
}

// main function
int main(int argc, const char** argv)
{
    

    // activate software
    mj_activate("mjkey.txt");

    
    // load and compile model
    char error[1000] = "Could not load binary model";

    // check command-line arguments
    if (argc < 2)
        m = mj_loadXML(filename, 0, error, 1000);

    else
        if (strlen(argv[1]) > 4 && !strcmp(argv[1] + strlen(argv[1]) - 4, ".mjb"))
            m = mj_loadModel(argv[1], 0);
        else
            m = mj_loadXML(argv[1], 0, error, 1000);
    if (!m)
        mju_error_s("Load model error: %s", error);

    // make data
    d = mj_makeData(m);


    // init GLFW
    if (!glfwInit())
        mju_error("Could not initialize GLFW");

    // create window, make OpenGL context current, request v-sync
    GLFWwindow* window = glfwCreateWindow(1244, 700, "Demo", NULL, NULL);
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    // initialize visualization data structures
    mjv_defaultCamera(&cam);
    mjv_defaultOption(&opt);
    mjv_defaultScene(&scn);
    mjr_defaultContext(&con);
    mjv_makeScene(m, &scn, 2000);                // space for 2000 objects
    mjr_makeContext(m, &con, mjFONTSCALE_150);   // model-specific context

    // install GLFW mouse and keyboard callbacks
    glfwSetKeyCallback(window, keyboard);
    glfwSetCursorPosCallback(window, mouse_move);
    glfwSetMouseButtonCallback(window, mouse_button);
    glfwSetScrollCallback(window, scroll);

    //double arr_view[] = {89.608063, -11.588379, 5, 0.000000, 0.000000, 0.000000};
    double arr_view[] = { 90, -5, 5, 0.012768, -0.000000, 1.254336 };
    cam.azimuth = arr_view[0];
    cam.elevation = arr_view[1];
    cam.distance = arr_view[2];
    cam.lookat[0] = arr_view[3];
    cam.lookat[1] = arr_view[4];
    cam.lookat[2] = arr_view[5];

    // qpos is dim xq x 1 = 7 x 1; 3 translations + 4 quaternions
    
    // custom controller
    mjcb_control = mycontroller; // 무한 반복되는 함수

    fid = fopen(datafile, "w");
    init_save_data();

// 초기 각도 입력
    d->qpos[0] = -pi/3;   // joint 1 - HFE

    // d -> qpos[1] = pi/2;
    //d->qpos[1] = 0.5;   // joint 2 - KFE
    theta[0] = d->qpos[0];
    


    state_init(m,d,&dPendulum,1);
    

    // use the first while condition if you want to simulate for a period.
    while (!glfwWindowShouldClose(window)) // 주기 : 0.0167
    {
    

        // advance interactive simulation for 1/60 sec
        //  Assuming MuJoCo can simulate faster than real-time, which it usually can,
        //  this loop will finish on time for the next frame to be rendered at 60 fps.
        //  Otherwise add a cpu timer and exit this loop when it is time to render.
        
        mjtNum simstart = d->time;
        

        // 여기가 loop time 0.0001인 부분 -> update하는 부분을 넣으면 되는 부분임
        while (d->time - simstart < 1.0 / 60.0) 
        {
            
            mj_step(m, d);
            // printf("err : %f, err_old : %f, err-err_old : %f\n",err,err_old, err-err_old);        
        }
        // printf("%f \n ", d -> time - simstart);
        if (d->time >= simend) {
            fclose(fid);
            break;
        }

        // get framebuffer viewport
        mjrRect viewport = { 0, 0, 0, 0 };
        glfwGetFramebufferSize(window, &viewport.width, &viewport.height);

        // update scene and render
        //opt.frame = mjFRAME_WORLD;
        //cam.lookat[0] = d->qpos[0];
        //cam.lookat[1] = 0;
        //cam.lookat[2] = 0;
        mjv_updateScene(m, d, &opt, NULL, &cam, mjCAT_ALL, &scn);
        mjr_render(viewport, &scn, &con);
        //printf("{%f, %f, %f, %f, %f, %f};\n",cam.azimuth,cam.elevation, cam.distance,cam.lookat[0],cam.lookat[1],cam.lookat[2]);

        // swap OpenGL buffers (blocking call due to v-sync)
        glfwSwapBuffers(window);

        // process pending GUI events, call GLFW callbacks
        glfwPollEvents();

    }

    // free visualization storage
    mjv_freeScene(&scn);
    mjr_freeContext(&con);

    // free MuJoCo model and data, deactivate
    mj_deleteData(d);
    mj_deleteModel(m);
    mj_deactivate();

    // terminate GLFW (crashes with Linux NVidia drivers)
    #if defined(__APPLE__) || defined(_WIN32)
        glfwTerminate();
    #endif

    return 1;
}
