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
// double ctrl_old = 0; 
// double ctrl = 0 ; 0
// double dense_M[4] = {0};
// double M[4] = {0};
double err_old=  0; double err = 0;  
double q= 0; 
// double q_old1 = 0;
// double q_old2 = 0;
double ref = -pi/2;
// double qd = 0; double qd_old = 0;
// double qdd[3] = {0};
double Ts = 0.0001; double cutoff = 250;
double ctrl_input;

void mycontroller(const mjModel* m, mjData* d)
{
// single pendulum
    q = d->qpos[0]; err = ref - q;
    
    //PID loop
    pid1.set_gain(15,0,2);
    
    // control torque
    ctrl_input = pid1.compute_PID(err,err_old, Ts, cutoff);
    // 1st order
    // mj_fullM(m, dense_M, d->qM);
    // double M[ndof][ndof] = { 0 };   // 2x2 inertia matrix -> 알아서 계산해줌.
    // M[0][0] = dense_M[0];
    // M[0][1] = dense_M[1];
    // M[1][0] = dense_M[2];
    // M[1][1] = dense_M[3];
    // printf("%f %f \n", M[0][0], M[0][1]);
    // printf("%f %f \n", M[1][0], M[1][1]);
    // printf("*********** \n");


    d->ctrl[0] = ctrl_input;
    
    // pid1.set_gain(10,3,0);
    // pid_result[0] = pid1.compute_PID(ref, dPendulum.q,0.0001, 100);
    // d->ctrl[0] = pid_result[0];

    // ref[1] = ref[0]; // update
    // dPendulum.q[1] = dPendulum.q[0];
    
    // d->qpos[0] = -pi/3;
    // d->qpos[1] = pi/3;
    // d->ctrl[0] +=pi/10 ;
    
    // double dense_M[ndof * ndof] = { 0 };   // 4x1 inertia matrix
    // mj_fullM(m, dense_M, d->qM);
    // double M[ndof][ndof] = { 0 };   // 2x2 inertia matrix -> 알아서 계산해줌.
    // M[0][0] = dense_M[0];
    // M[0][1] = dense_M[1];
    // M[1][0] = dense_M[2];
    // M[1][1] = dense_M[3];
    // // printf("%f %f \n", M[0][0], M[0][1]);
    // // printf("%f %f \n", M[1][0], M[1][1]);
    // // printf("*********** \n");
    // pid1.set_gain(0.1,0.1,0);
    // pid2.set_gain(0.1,0.1,0);
    // update_state(m,d,&dPendulum);
    // // double r, double y, double dt, double Wc

    // printf("err : %f, err_old : %f, err-err_old : %f, dterm : %f\n",err,err_old, err-err_old, dterm);
    // pterm = kp*err;
    // dterm = kd*(2*(err-err_old)-(dt-2*t_const)*dterm_old)/(2*t_const+dt);
    // d->ctrl[0] = pterm + dterm;
    // err_old = err; dterm_old = dterm;
    

    // pid_result[0] = pid1.compute_PID(ref[0], dPendulum.q[0],0.0001, 100);
    // pid_result[1] = pid2.compute_PID(ref[1], dPendulum.q[1],0.0001, 100);
    // dPendulum.ctrl_input[0] = pid_result[0];
    // dPendulum.ctrl_input[1] = pid_result[1];
    // mju_mulMatVec(torque_input,dense_M,dPendulum.ctrl_input,2,2);
    
    // //dob 설계
    // double d_hat[ndof] = { 0};
    // double m_n[ndof] = {0};

    // mju_mulMatVec(m_n, dense_M, dPendulum.qdd,2,2);
    // // mju_scl(m_n, m_n, -1, 2);
    // mju_sub(d_hat, torque_input, m_n, 2);   
    // mju_sub(torque_input,torque_input, d_hat,2);
    
    // d->ctrl[0] = torque_input[0];
    // d->ctrl[1] = torque_input[1];
    // printf("torque1 : %f,torque1 : %f\n", torque_input[0], torque_input[1]);
    // printf("error1 : %f, error2 = %f\n", ref[0]- d->qpos[0], ref[1] -d->qpos[1]);




    
    // double qddot[ndof] = { 0 };
    // qddot[0] = d->qacc[0];
    // qddot[1] = d->qacc[1];

    // double h[ndof] = { 0 };
    // h[0] = d->qfrc_bias[0];
    // h[1] = d->qfrc_bias[1];

    // double lhs[ndof] = { 0 };
    // mju_mulMatVec(lhs, dense_M, qddot, 2, 2);   // lhs = M*qddot
    // // lhs[0] = lhs[0] + h[0]; // lhs = M*qddot + h
    // // lhs[1] = lhs[1] + h[1];
    // // d->qfrc_applied[0] = lhs[0];
    // // d->qfrc_applied[1] = lhs[1];
    // d->ctrl[0] = lhs[0];
    // d->ctrl[1] = lhs[1];
    


    // PID pid_1;
    // PID pid_2;
    // pid_1.set_gain(17,4,0);
    // pid_2.set_gain(7,1,0);
    // // for(int i= 0; i<2; i++)
    // // {
    // ctrl_input[0] = pid_1.compute_PID(ref[0], d->qpos[0],0.0001, 100);
    // d->ctrl[0] = ctrl_input[0]; 
    // ctrl_input[1] = pid_2.compute_PID(ref[1], d->qpos[1],0.0001, 0.1);
    // d->ctrl[1] = ctrl_input[1]; 
    // mj_fullM()


    // }

    if (loop_index % data_frequency == 0) {     // loop_index를 data_frequency로 나눈 나머지가 0이면 데이터를 저장.
        save_data(m, d);
        // ref[0] += pi/100;
        // ref[1] += pi/200;
        // printf("error 1 : %f,2 : %f \n", ref[0] - d->qpos[0], ref[1] - d->qpos[1]);

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

    fid = fopen(datapath, "w");
    init_save_data();

// 초기 각도 입력
    d->qpos[0] = pi/3;   // joint 1 - HFE
    // d -> qpos[1] = pi/2;
    //d->qpos[1] = 0.5;   // joint 2 - KFE
    


    state_init(m,d,&dPendulum,1);

    ;

    // use the first while condition if you want to simulate for a period.
    while (!glfwWindowShouldClose(window))
    {
        // advance interactive simulation for 1/60 sec
        //  Assuming MuJoCo can simulate faster than real-time, which it usually can,
        //  this loop will finish on time for the next frame to be rendered at 60 fps.
        //  Otherwise add a cpu timer and exit this loop when it is time to render.
        mjtNum simstart = d->time;
        while (d->time - simstart < 1.0 / 60.0)
        {
            mj_step(m, d);
            printf("err : %f, err_old : %f, err-err_old : %f\n",err,err_old, err-err_old);


            err_old = err;
            pid1.update_PID();


            // Drag Force = -c*v^2*unit_vector(v); v = sqrt(vx^2+vy^2+vz^2);
            // vector(v) = vx i + vy j + vz k;
            // unit vector = vector(v)/v;
            // fx = -c*v*vx; fy = -c*v*vy; fz = -c*v*vz;

            //double vx, vy, vz;
            //vx = d->qvel[0]; vy = d->qvel[1]; vz = d->qvel[2];

            //double v;
            //v = sqrt(vx^2 + vy^2 + vz^2); <-- 이렇게 하면 안되는듯 
            //v = sqrt(vx * vx + vy * vy + vz * vz);

            //double fx, fy, fz, c;
            //c = 0.5;
            //fx = -c * v * vx;
            //fy = -c * v * vy;
            //fz = -c * v * vz;

            //d->qfrc_applied[0] = fx;
            //d->qfrc_applied[1] = fy;
            //d->qfrc_applied[2] = fz;
            
        }
        
        // if (d->time >= simend) {
        //     fclose(fid);
        //     break;
        // }

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
