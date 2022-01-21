#include <bits/stdc++.h>
#include "OP_new.h"
#include "programs_new.h"
#include "sun_model_new.h"
#include "ephemerides_new.h"
#include "frame_conversion_new.h"
#include "estimation.h"
#include "disturbance_torque.h"
#include "bdot.h"
#include "integrator.h"
#include "pid.h"
#define pi 3.14159265359

using namespace std;
ifstream sim_file("sim_state.txt");
ifstream sunbody_sil("sun_body2.txt");
ifstream magbody_sil("mag_body2.txt");
ifstream sunorb_sil("sun_orb2.txt");
ifstream magorb_sil("mag_orb2.txt");
ifstream mag_ref("mag_out.txt");
ofstream magfield_eci("mag_eci.txt");
ofstream magfield_orbit("mag_orbit.txt");
ofstream magfield_body("mag_body.txt");
ofstream sunvec_eci("sun_eci.txt");
ofstream sunvec_body("sun_body.txt");
ofstream sunvec_orbit("sun_orbit.txt");
ofstream lla_pos("../py-igrf/release/lla_pos.txt");
ofstream wbi("w_bi.txt");
ofstream wri("w_ri.txt");
ofstream w_err("w_error.txt");
ofstream lla("lla_position.txt");
ofstream q_err("q_error.txt");
ofstream q_err_1("q_error_1.txt");
ofstream q_bi("q_bi.txt");
ofstream q_bi_1("q_bi_1.txt");
ofstream q_ri("q_ri.txt");
ofstream q_br("q_br.txt");
ofstream tc("control_torque.txt");

double **adcs_engine(int time,double**moi_satellite, double**prev_req_var){
    
    double **position, **velocity, **state_sim, **position_ecef, **position_lla, **state_onboard, **state_onboard_new, **keplerian_elements; 
    double initial_julian_date, current_julian_date, gstangle;
    double **sun_out, **sun_eci, **sun_orbit, **sun_body, **magnetic_field_py, **mag_body, **mag_eci_frame, **mag_orbit;
    double **trans, **moi_sat_inverse;
    double **w_ri, **w_error, **w_bi, **w_bi_updated;
    double **quaternion_br, **quaternion_bi, **quaternion_bi_1, **quaternion_bi_1_updated,**quaternion_ri, **quaternion_ri_prev, **quaternion_error, **quaternion_error_1, **quaternion_error_integral;
    double **control_torque, **disturbance_torques;
    double** sun_orbit_sil, **sun_body_sil, **mag_orbit_sil, **mag_body_sil;

    int i,j;
    double step = 1.0;

    if(time%1000 ==0)
    cout<<time<<endl;
    
    position = getzeromatrix(3,1);
    velocity = getzeromatrix(3,1);
    state_onboard = getzeromatrix(6,1);
    magnetic_field_py = getzeromatrix(3,1);
    state_sim = getzeromatrix(6,1);
    sun_eci = getzeromatrix(3,1);
    quaternion_bi_1 = getzeromatrix(4,1);
    quaternion_error_integral = getzeromatrix(4,1);
    quaternion_ri_prev = getzeromatrix(4,1);
    sun_body_sil = getzeromatrix(3,1);
    sun_orbit_sil = getzeromatrix(3,1);
    mag_body_sil = getzeromatrix(3,1);
    mag_orbit_sil = getzeromatrix(3,1);
    w_bi = getzeromatrix(3,1);
    
    for(i=0;i<4;i++){
        quaternion_bi_1[i][0] = prev_req_var[i][0];
        quaternion_ri_prev[i][0] = prev_req_var[i+4][0];
        quaternion_error_integral[i][0] = prev_req_var[i+8][0];
        
    }

    for(i=0;i<3;i++){
        w_bi[i][0] = prev_req_var[18+i][0];
    }

    string temp;
    getline(sim_file,temp);
    getline(mag_ref,temp);
    getline(magbody_sil,temp);
    getline(magorb_sil,temp);
    getline(sunbody_sil,temp);
    getline(sunorb_sil,temp);
    
    for(i=0;i<6;i++){
        sim_file >> state_sim[i][0];
        state_onboard[i][0] = prev_req_var[i+12][0];
    }

    for(i=0;i<3;i++){
        sunorb_sil >> sun_orbit_sil[i][0];
        sunbody_sil >> sun_body_sil[i][0];
        magorb_sil >> mag_orbit_sil[i][0];
        magbody_sil >> mag_body_sil[i][0];
    }

    state_onboard_new = orbit_propagate(state_onboard,time-1,time,0.01);
    keplerian_elements = cartesian_to_keplerian(state_onboard_new);

    for(i=0;i<3;i++){
        position[i][0] = state_onboard_new[i][0];
        velocity[i][0] = state_onboard_new[i+3][0];
    }

    initial_julian_date = getJulianDate(2014,3,17,23,25,30);
    current_julian_date = initial_julian_date + (time/86400.0);
    gstangle =  jdut2gmst(current_julian_date);
    
    position_ecef = eci_to_ecef(position,gstangle);

    position_lla = ecef_to_lla(position_ecef);

    trans = eci_orbit(keplerian_elements);
    
    sun_out = sun_model(current_julian_date,state_onboard_new);
   
    for(i=0;i<3;i++){
        sun_eci[i][0] = sun_out[i+1][0];
    }

    sun_orbit = matrixmultiply(trans,sun_eci,3,3,1);
    sun_body = eci_body(sun_eci, quaternion_bi_1);

    for(i=0;i<3;i++){
        mag_ref >> magnetic_field_py[i][0];
    }

    double **magnetic_field_py_tesla = MatrixScalarMultiply(magnetic_field_py,1e-9,3,1);
    mag_eci_frame = ned_to_ecef_to_eci(magnetic_field_py_tesla,position_lla,gstangle);
    mag_orbit = matrixmultiply(trans,mag_eci_frame,3,3,1); 
    mag_body = eci_body(mag_eci_frame, quaternion_bi_1);

    double **unit_sun_body, **unit_mag_body, **unit_sun_orbit, **unit_mag_orbit, **quaternion_rb;

    unit_mag_body = Unit_vector(mag_body,3,1);
    unit_mag_orbit = Unit_vector(mag_orbit,3,1);
    unit_sun_body = Unit_vector(sun_body,3,1);
    unit_sun_orbit = Unit_vector(sun_orbit,3,1);
   
    quaternion_rb = quest(0.5,0.5,unit_sun_orbit,unit_sun_body,unit_mag_orbit,unit_mag_body); 
    quaternion_br = quatinv(quaternion_rb);
    quaternion_ri = compute_q_ri(trans,quaternion_ri_prev);
    quaternion_bi = quatmultiply(quaternion_br, quaternion_ri);

    double **q_bi_inv, **q_bi_1_inv;
    q_bi_inv = quatinv(quaternion_bi);
    q_bi_1_inv = quatinv(quaternion_bi_1);
    quaternion_error = quatmultiply(q_bi_inv,quaternion_ri);
    quaternion_error_1 = quatmultiply(q_bi_1_inv,quaternion_ri);

    w_ri = step_differentiate(step, quaternion_ri, quaternion_ri_prev);
    w_error = MatrixSubtract(w_ri, w_bi, 3, 1);
    double w_bi_norm = norm(w_bi);

    if (w_bi_norm > 0.0872665){
        control_torque = detumble(mag_body,w_bi,moi_satellite); 
    }
    else{
        control_torque = tc_pid(moi_satellite, quaternion_error, w_error, quaternion_error_integral);
    }

    disturbance_torques = calc_disturbance_torque(keplerian_elements, state_onboard_new, quaternion_br, sun_body);
    // disturbance_torques = getzeromatrix(3,1);

    moi_sat_inverse = InverseofMatrix(moi_satellite,3,3);

    w_bi_updated = integrate_omega(moi_satellite,moi_sat_inverse,control_torque,disturbance_torques,w_bi, 1, 0.2); 
    quaternion_bi_1_updated = integrate_quaternion(quaternion_bi_1, w_bi_updated, 1, 0.2);

    //storing previously required variable in an array
    for(i=0;i<4;i++){
        prev_req_var[i][0] = quaternion_bi_1_updated[i][0];
        prev_req_var[i+4][0] = quaternion_ri[i][0];
        prev_req_var[i+8][0] = quaternion_error_integral[i][0];
    }
    for(i=0;i<6;i++){
        prev_req_var[i+12][0] = state_onboard_new[i][0];
    }
    for(i=0;i<3;i++){
        prev_req_var[i+18][0] = w_bi_updated[i][0];
    }
    //printing to file
    lla<<time<<" "<<position_lla[0][0]<<" "<<position_lla[1][0]<<" "<<position_lla[2][0]<<endl;
    lla_pos<<position_lla[0][0]<<" "<<position_lla[1][0]<<" "<<position_lla[2][0]<<endl;
    sunvec_eci<<time<<" "<<sun_eci[0][0]<<" "<<sun_eci[1][0]<<" "<<sun_eci[2][0]<<endl;
    sunvec_orbit<<time<<" "<<sun_orbit[0][0]<<" "<<sun_orbit[1][0]<<" "<<sun_orbit[2][0]<<endl;
    sunvec_body<<time<<" "<<sun_body[0][0]<<" "<<sun_body[1][0]<<" "<<sun_body[2][0]<<endl;
    magfield_eci<<time<<" "<<mag_eci_frame[0][0]<<" "<<mag_eci_frame[1][0]<<" "<<mag_eci_frame[2][0]<<endl;
    magfield_body<<time<<" "<<mag_body[0][0]<<" "<<mag_body[1][0]<<" "<<mag_body[2][0]<<endl;
    magfield_orbit<<time<<" "<<mag_orbit[0][0]<<" "<<mag_orbit[1][0]<<" "<<mag_orbit[2][0]<<endl;
    q_ri<<time<<" "<<quaternion_ri[0][0]<<" "<<quaternion_ri[1][0]<<" "<<quaternion_ri[2][0]<<" "<<quaternion_ri[3][0]<<endl;
    q_br<<time<<" "<<quaternion_br[0][0]<<" "<<quaternion_br[1][0]<<" "<<quaternion_br[2][0]<<" "<<quaternion_br[3][0]<<endl;
    q_bi<<time<<" "<<quaternion_bi[0][0]<<" "<<quaternion_bi[1][0]<<" "<<quaternion_bi[2][0]<<" "<<quaternion_bi[3][0]<<endl;
    q_bi_1<<time<<" "<<quaternion_bi_1[0][0]<<" "<<quaternion_bi_1[1][0]<<" "<<quaternion_bi_1[2][0]<<" "<<quaternion_bi_1[3][0]<<endl;
    q_err<<time<<" "<<quaternion_error[0][0]<<" "<<quaternion_error[1][0]<<" "<<quaternion_error[2][0]<<" "<<quaternion_error[3][0]<<endl;
    q_err_1<<time<<" "<<quaternion_error_1[0][0]<<" "<<quaternion_error_1[1][0]<<" "<<quaternion_error_1[2][0]<<" "<<quaternion_error_1[3][0]<<endl;
    tc<<time<<" "<<control_torque[0][0]<<" "<<control_torque[1][0]<<" "<<control_torque[2][0]<<endl;
    wbi<<time<<" "<<w_bi[0][0]*(180.0/pi)<<" "<<w_bi[1][0]*(180.0/pi)<<" "<<w_bi[2][0]*(180.0/pi)<<endl;
    wri<<time<<" "<<w_ri[0][0]*(180.0/pi)<<" "<<w_ri[1][0]*(180.0/pi)<<" "<<w_ri[2][0]*(180.0/pi)<<endl;
    w_err<<time<<" "<<w_error[0][0]*(180.0/pi)<<" "<<w_error[1][0]*(180.0/pi)<<" "<<w_error[2][0]*(180.0/pi)<<endl;

    free_variable(position,3);
    free_variable(velocity,3);
    free_variable(state_sim,6);

    free_variable(magnetic_field_py,3);
    free_variable(magnetic_field_py_tesla,3);

    free_variable(quaternion_bi,4);
    free_variable(quaternion_rb,4);
    free_variable(quaternion_br,4);
    free_variable(quaternion_ri,4);
    free_variable(quaternion_bi_1_updated,4);
    free_variable(quaternion_bi_1,4);
    free_variable(quaternion_ri_prev,4);

    free_variable(quaternion_error,4);
    free_variable(quaternion_error_1,4);
    free_variable(quaternion_error_integral,4);

    free_variable(position_lla,3);
    free_variable(position_ecef,3);

    free_variable(sun_out,4);
    free_variable(sun_eci,3);
    free_variable(sun_orbit,3);
    free_variable(sun_body,3);
    free_variable(mag_body,3);

    free_variable(unit_sun_orbit,3);
    free_variable(unit_sun_body,3);
    free_variable(unit_mag_body,3);
    free_variable(mag_eci_frame,3);
    free_variable(mag_orbit,3);
    free_variable(unit_mag_orbit,3);

    free_variable(control_torque,3);
    free_variable(disturbance_torques,3);

    free_variable(state_onboard,6);
    free_variable(keplerian_elements,6);

    free_variable(w_bi,3);
    free_variable(w_error,3);
    free_variable(w_ri,3);

    free_variable(q_bi_1_inv,4);
    free_variable(q_bi_inv,4);

    free_variable(moi_sat_inverse,3);
    free_variable(trans,3);

    free_variable(sun_body_sil,3);
    free_variable(sun_orbit_sil,3);
    free_variable(mag_body_sil,3);
    free_variable(mag_orbit_sil,3);

    return prev_req_var;
}

int main(){

    double **state_onboard,**keplerian_elements,**orientation,**w_bi,**moi_satellite,**q_ri_init,**q_bi_1,**prev_required_var, **trans;

    state_onboard = getzeromatrix(6,1);
    orientation = getzeromatrix(3,1);
    w_bi = getzeromatrix(3,1);
    moi_satellite = getzeromatrix(3,3);
    prev_required_var = getzeromatrix(21,1);
    q_bi_1 = getzeromatrix(4,1);

    // state_onboard[0][0] =  -5236.84633;
    // state_onboard[1][0] =  4124.17773;
    // state_onboard[2][0] =  -1262.94137;

    // state_onboard[3][0] =  -3.86204515;
    // state_onboard[4][0] =  -3.12048032;
    // state_onboard[5][0] =  5.83839029;

    state_onboard[0][0] = 1238.488583514314;
    state_onboard[1][0] = 7080.353211119648;
    state_onboard[2][0] = 14.586295240982;

    state_onboard[3][0] = 1.111578762356251;
    state_onboard[4][0] = -0.204702311583077;
    state_onboard[5][0] = 7.364195645361482;

    keplerian_elements = cartesian_to_keplerian(state_onboard);

    // orientation[0][0] = 0.052359;
    // orientation[1][0] = 0.069813;
    // orientation[2][0] = 0.087266;

    orientation[0][0] = -1.5842;
    orientation[1][0] = 1.4185;
    orientation[2][0] = 0.1864;

    trans = eci_orbit(keplerian_elements);

    double **quat_trans;

    quat_trans = trans_to_quaternion(trans);

    q_ri_init  = quatinv(quat_trans);

    for (int i=0; i<4; i++){
        q_bi_1[i][0] = q_ri_init[i][0];
    }

    // cout<<q_br_init[0][0]<<" "<<q_br_init[1][0]<<" "<<q_br_init[2][0]<<" "<<q_br_init[3][0]<<endl;

    for(int i=0;i<4;i++){
        prev_required_var[i][0] = q_bi_1[i][0];
        prev_required_var[i+4][0] = q_ri_init[i][0];
    }
    for(int i=0;i<6;i++){
        prev_required_var[i+12][0] = state_onboard[i][0];
    }


    // -0.0572483 0.00999538 0.0089191
    w_bi[0][0] = -0.0572483;
    w_bi[1][0] = 0.00999538;
    w_bi[2][0] = 0.0089191;

    moi_satellite[0][0] = 4.46e-3;
    moi_satellite[1][1] = 1.16e-2;
    moi_satellite[2][2] = 1.15e-2;

    for(int i=0;i<3;i++){
        prev_required_var[i+18][0] = w_bi[i][0];
    }

    int t = 3;
    for(int time=1;time<t;time++) {
       prev_required_var = adcs_engine(time,moi_satellite,prev_required_var);
    }

    free_variable(state_onboard,6);
    free_variable(keplerian_elements,6);
    free_variable(prev_required_var,21);
    free_variable(q_bi_1,4);
    free_variable(orientation,3);
    free_variable(moi_satellite,3);
    free_variable(w_bi,3);
    free_variable(q_ri_init,4);
    free_variable(quat_trans,4);
    free_variable(trans,3);

    return 0;
}