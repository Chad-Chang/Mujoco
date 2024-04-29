// #pragma once

#include "funcs.h" // 사용자 정의 헤더파일 <>하면 안됨.

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
int vel_count = 0;
int acc_count = 0;
double alpha[3] = {0};
double beta[3] = {0};
void mycontroller(const mjModel* m, mjData* d)  // 제어주기 0.000025임
{
// single pendulum
    if(loop_index % contol_loop ==0) // sampling time 0.0001
    {   
        vel_count += 1;
        acc_count += 1;
        perturb = amplitude_perturb*sin(dist_freq*2*pi*d->time);
        // loop_tcheck();
        err = ref - theta[0];
        pid1.set_gain(7,0,1.3);
            // control torque 
        u_c = pid1.compute_PID(err,err_old, Ts, cutoff);

        
        u_d[0] = u_c - d_hat[0];
        

        d->ctrl[0] = u_d[0] + perturb;
        err_old = err;
        pid1.update_PID();


        theta[0] = d->qpos[0];

        theta_vel[0] = tustin_derivative(theta[0],theta[1], theta_vel[1], LPF_freq);

        theta_acc[0] = tustin_derivative(theta_vel[0],theta_vel[1], theta_acc[1], LPF_freq);

        alpha[0] = lowpassfilter(u_d[0], u_d[1], alpha[1], LPF_freq);
        beta[0] = lowpassfilter(alpha[0],alpha[1], beta[1], LPF_freq);
        
        
        
        // if(vel_count>2)
        // {
        
        // theta_vel[0] = (tau*(theta[0]-theta[1])-(Ts-2*tau)*theta_vel[1])/(2*tau+Ts);
        // }

        // if(acc_count>3)
        // {
        
        // theta_acc[0] = (tau*(theta_vel[0]-theta_vel[1])-(Ts-2*tau)*theta_acc[1])/(2*tau+Ts);
        // }
        // theta_acc[0] = (4*(theta[2]-2*theta[1]+theta[0]) + (8*pow(tau,2)-2*pow(Ts,2))*theta[1]-(pow(Ts,2)-4*tau*Ts+4*pow(tau,2)))/(pow(Ts,2)+4*tau*Ts+4*pow(tau,2));
        // theta_acc[0] = 

        
    // 미분마다 low pass filter 씌운 놈
        // theta_vel[2] = theta_vel[1];
        // theta_vel[1] = theta_vel[0];
        // theta_vel[0] = (2*(theta[0]-theta[1])-(Ts-2*tau)*theta_vel[1])/(2*tau + Ts);

        // theta_acc[2] = theta_acc[1];
        // theta_acc[1] = theta_acc[0];
        // theta_acc[0] = (-(Ts-2*tau)*theta_acc[1]+2*(theta_vel[0]-theta_vel[1]))/(2*tau +Ts);

    // 필터랑 미분 한번에 한놈

    // 그냥 미분만 해준놈(필터 없음)
        // if(vel_count > 2)
        // {
        //     theta_vel[2] = theta_vel[1]; theta_vel[1] = theta_vel[0];
        //     theta_vel[0] = (2/Ts)*(theta[0]-theta[1]) - theta_vel[1];
        // }
        
        // if(acc_count > 3)
        // {
        //     theta_acc[2] = theta_acc[1]; 
        //     theta_acc[1] = theta_acc[0];
        //     theta_acc[0] = 4/pow(Ts,2)*(theta[0]-2*theta[1]+theta[2])-(2/Ts)*(theta_vel[1]-theta_vel[2])-theta_acc[1];
        // }
        // if(acc_count < 3)
        //     printf("before theta_acc = %f\n",theta_acc[0]);
        // printf("dist_error  = %f\n", perturb - d_hat[0]);
        // printf("IO error = %f\n",ref-d->qpos[0]);



        mj_fullM(m, dense_M, d->qM);
        double M[ndof][ndof] = { 0 };   // 2x2 inertia matrix -> 알아서 계산해줌.
        M[0][0] = dense_M[0];
        M[0][1] = dense_M[1];
        M[1][0] = dense_M[2];
        M[1][1] = dense_M[3];


        
        // if(acc_count >= 3)
        //     printf("After theta_acc = %f\n",theta_acc[0]);

    
        // printf("M = %f\n", M[0][0]);

        // d_hat[1]= d_hat[0];
        if(dob_switch)
                {
                    d_hat[0] = M[0][0]*theta_acc[0]-beta[0];
                    
                    theta[2] = theta[1]; theta[1] = theta[0];
                    theta_vel[2] = theta_vel[1]; theta_vel[1] = theta_vel[0];
                    theta_acc[2] = theta_acc[1]; theta_acc[1] = theta_acc[0];
                    alpha[1] = alpha[0];
                    beta[1] = beta[0];
                    u_d[2] = u_d[1]; u_d[1] = u_d[0];
                    
                }
            // else
            //     printf("d_ = %f\n", d_hat[0]);

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
    d->qpos[0] = pi/3;   // joint 1 - HFE
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
