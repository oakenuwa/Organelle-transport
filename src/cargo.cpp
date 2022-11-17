/*------------------------------------------------------------------
 cargo.cpp : object describing a cargo-motor complex
 
 Copyright (C) 2016
 Created by: Simon Freedman, Shiladitya Banerjee, Glen Hocky, Aaron Dinner
 Contact: dinner@uchicago.edu

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version. See ../LICENSE for details.
-------------------------------------------------------------------*/

#include "globals.h"
#include "motor.h"
#include "cargo.h"
#include "filament_ensemble.h"
//#include "actin.h"

//motor class

cargo::cargo(array<double, 2> myfov, double delta_t, double temp,
            filament_ensemble * network,
             double vis, string BC)
{
    //vector<motor *> n_cargo;
    fov = myfov;
    gamma = 0;
    tMove=0;//10;
    f_network=network;
    bc = BC;
    
    ke = 0;
    pe = 0;
    cdt = delta_t;
    cdamp = (6*pi*vis*0.8); //*mld
    ctemperature = temp;
    cbd_prefactor = sqrt(ctemperature/(2*cdamp*cdt));
    c_prv_rnd_x = 0;
    c_prv_rnd_y = 0;
}

cargo::cargo(int mpercargo, array<double, 2> myfov, double delta_t, double temp,
    double mlen, filament_ensemble * network, double v0, double stiffness, double max_ext_ratio,
    double ron, double roff, double rend,
    double fstall, double rcut,
    double vis, vector<array<double,3> > positions, string BC) {

    fov = myfov;
    mld = mlen;
    gamma = 0;
    tMove=0;//10;
    f_network=network;
    v = v0;
    bc = BC;
    
    ke = 0;
    pe = 0;
    cdt = delta_t;
    cdamp = (6*pi*vis*0.8); //*mld
    ctemperature = temp;
    cbd_prefactor = sqrt(ctemperature/(2*cdamp*cdt));
    int motor_indicator = 1;
    
    
    int nm = mpercargo;
    //cout<<"\nDEBUG: Number of motors per cargo:"<<nm<<"\n";

    double alpha = 1, motorx, motory, mang;
    array<double, 3> motor_pos;
    cmotorx = rng(-0.5*(fov[0]*alpha-(2*mld)),0.5*(fov[0]*alpha-(2*mld)));
    cmotory = rng(-0.5*(fov[1]*alpha-(2*mld)),0.5*(fov[1]*alpha-(2*mld)));
    
    c_prv_rnd_x = 0;
    c_prv_rnd_y = 0;
    for (int i=0; i< nm; i++) {

    if ((unsigned int)i < positions.size()){
        motorx = positions[i][0];
        motory = positions[i][1];
        mang   = positions[i][2];
    }else{
        //motorx = rng(-0.5*(fov[0]*alpha-mld),0.5*(fov[0]*alpha-mld));
        //motory = rng(-0.5*(fov[1]*alpha-mld),0.5*(fov[1]*alpha-mld));
        mang   = rng(0,2*pi);
        cout<<"my angle is: "<<mang<<endl;
    }
    motor_pos = {cmotorx, cmotory, mang};
    n_cargo.push_back(new motor(motor_pos, mld, f_network,{0, 0}, {-1,-1}, {-1,-1}, fov, delta_t, temp,
                v0, stiffness, max_ext_ratio, ron, roff, rend, fstall, rcut, vis, BC, motor_indicator));
    }
    
}

//create a function that creates a cargo object from a list of cargo positions
//create a vector of cargos, then add
void cargo::add_motor(motor * m){
    //cout<<"Motor m added"<<endl;
    n_cargo.push_back(m);
    array<double,2> p = m->get_hx();
    array<double,2> v = m->get_hy();
    cmotorx = p[0];
    cmotory = v[0];
}
/*cargo::cargo( array<double, 4> pos,
        double mlen, filament_ensemble * network,
        array<int, 2> mystate,
        array<int, 2> myfindex,
        array<int, 2> mylindex,
        array<double, 2> myfov,
        double mycargoindex,
        double delta_t,
        double temp,
        double v0,
        double stiffness,
        double max_ext_ratio,
        double ron, double roff, double rend,
        double fstall, double rcut,
        double vis, string bc) {
    
    vs          = v0;
    mk          = stiffness;
    
    stall_force = fstall;
    temperature = temp;

    max_bind_dist = rcut;
    
    mld         = mlen;
    dt          = delta_t;
    kon         = ron*dt;
    koff        = roff*dt;
    kend        = rend*dt;
    
    state       = mystate;
    f_index     = myfindex; //filament index for each head
    l_index     = mylindex; //link index for each head
    fov         = myfov;
    BC          = bc;
    actin_network = network;
    damp        =(6*pi*vis*mld);
    bd_prefactor= sqrt(temperature/(2*damp*dt));
    
    ********for FENE springs*********
    max_ext     = max_ext_ratio*mlen;
    eps_ext     = 0.01*max_ext;
    ********************************

    shear       = 0;
    tension     = 0;
    force       = {0,0}; // force on the spring
    kinetic_energy = 0;
    pos_a_end = {0, 0}; // pos_a_end = distance from pointy end -- by default 0
                        // i.e., if l_index[hd] = j, then pos_a_end[hd] is the distance to the "j+1"th actin

    
    array<double, 2> posH0 = boundary_check(0, pos[0], pos[1]);
    array<double, 2> posH1 = boundary_check(1, pos[0]+pos[2], pos[1]+pos[3]);
    hx[0] = posH0[0];
    hy[0] = posH0[1];
    hx[1] = posH1[0];
    hy[1] = posH1[1];
   
    //force can be non-zero and angle is determined from disp vector
    this->update_angle();
    this->update_force();
    
    ldir_bind[0] = {0,0};
    ldir_bind[1] = {0,0};
    bind_disp[0] = {0,0};
    bind_disp[1] = {0,0};

    at_barbed_end = {false, false};

    if (state[0] == 1){
        pos_a_end[0] = dist_bc(BC, actin_network->get_end(f_index[0], l_index[0])[0] - hx[0],
                                   actin_network->get_end(f_index[0], l_index[0])[1] - hy[0], fov[0], fov[1], 0);
        ldir_bind[0] = actin_network->get_direction(f_index[0], l_index[0]);

    }
    if (state[1] == 1){
        pos_a_end[1] = dist_bc(BC, actin_network->get_end(f_index[1], l_index[1])[0] - hx[1],
                                   actin_network->get_end(f_index[1], l_index[1])[1] - hy[1], fov[0], fov[1], 0);
        ldir_bind[1] = actin_network->get_direction(f_index[1], l_index[1]);
    }

    prv_rnd_x = {0,0};
    prv_rnd_y = {0,0};

}
*/

 cargo::~cargo(){};

//return motor state with a given head number

array<int, 2> cargo::get_states(int motor)
{
    return n_cargo[motor]->get_states();
}

array<double, 2> cargo::get_hx(int motor)
{
    return n_cargo[motor]->get_hx();
}


array<double, 2> cargo::get_hy(int motor)
{
    return n_cargo[motor]->get_hy();
}


string cargo::get_BC(int motor)
{
    return n_cargo[motor]->get_BC();
}


//metropolis algorithm with rate constant
//NOTE: while fl_idx doesn't matter for this xlink implementation, it does for "spacers"



 
 void cargo::motor_walk(double t)
 {
     //cout<<"\ncargo_size: "<<n_cargo.size()<<endl;
     this->check_broken_filaments();
     int ncargo_sz = int(n_cargo.size());
     //bool attached;
     //#pragma omp parallel for
     double totforcex = 0;
     double totforcey = 0;
     for (int j=0; j<ncargo_sz; j++)
    {
        cout<<"motor "<<j<<endl;
       double forcex = n_cargo[j]->get_force()[0];
       double forcey = n_cargo[j]->get_force()[1];
        totforcex += forcex;
        totforcey += forcey;
    }
     //cout<<"force x: "<<totforcex<<" force y: "<<totforcey<<endl;
     double new_rnd_cx= rng_n(0,1), new_rnd_cy = rng_n(0,1);
     //cout<<"new rnd x: "<<new_rnd_cx<<" new rnd y: "<<new_rnd_cy<<endl;
     double vx =  totforcex/cdamp + cbd_prefactor*(new_rnd_cx + c_prv_rnd_x); //pow(-1,0)*
     double vy =  totforcey/cdamp + cbd_prefactor*(new_rnd_cy + c_prv_rnd_y); //pow(-1,0)*
     //cout<<"\nvx: "<<vx<<" vy: "<<vy<<endl;
     //kinetic_energy = vx*vx + vy*vy;
     array<double, 2> cpos = cargo_boundary_check(cmotorx + vx*cdt, cmotory + vy*cdt);
     cmotorx = cpos[0];
     cmotory = cpos[1];
    // cout<<"\npos x: "<<cpos[0]<<" pos y: "<<cpos[1]<<endl;
     for (int j = 0; j<ncargo_sz; j++)
     {
         n_cargo[j]->set_hx(0,cmotorx);
         n_cargo[j]->set_hy(0,cmotory);
         
     }
     
     c_prv_rnd_x = new_rnd_cx;
     c_prv_rnd_y = new_rnd_cy;
     
     for (int i=0; i<ncargo_sz; i++) {
         bool attached;
     //    if(i==0) cout<<"\nDEBUG: motor_walk: using "<<omp_get_num_threads()<<" cores";
         array<int, 2> s = n_cargo[i]->get_states();
         
         if (t >= tMove){
             if (s[1] == 1)
                 n_cargo[i]->step_onehead_motors(1);
             else if (s[1] == 0)
             {
                 //cout<<"head one about to attach"<<endl;
                 attached = n_cargo[i]->attach_motors(1);
                 if (!attached)
                     n_cargo[i]->brownian_relax(1);
             }

             n_cargo[i]->update_angle();
             n_cargo[i]->update_force();
             //n_motors[i]->update_force_fraenkel_fene();
             n_cargo[i]->actin_update();
         }
     
     }
     
     cout<<"cargo walk done at time "<<t<<" secs"<<endl;
    // this->update_energies();
     
 } 

array<double, 2> cargo::cargo_boundary_check(double x, double y)
{
    return pos_bc(bc, f_network->get_delrx(), cdt, fov, {(x - cmotorx)/cdt, (y - cmotory)/cdt}, {x, y});
}

void cargo::check_broken_filaments()
{
    vector<int> broken_filaments = f_network->get_broken();
    array<int, 2> f_index;
    
    for (unsigned int i = 0; i < broken_filaments.size(); i++){
        
        for(unsigned int j = 0; j < n_cargo.size(); j++){
            
            f_index = n_cargo[j]->get_f_index();

            if(f_index[0] == broken_filaments[i]){
                n_cargo[j]->detach_head_without_moving(0);
                //cout<<"\nDEBUG: detaching head 0 of motor "<<j<<" from broken filament "<<broken_filaments[i];
            }

            if(f_index[1] == broken_filaments[i]){
                n_cargo[j]->detach_head_without_moving(1);
                //cout<<"\nDEBUG: detaching head 1 of motor "<<j<<" from broken filament "<<broken_filaments[i];
            }
        }
    }


}


void cargo::kill_heads(int hd)
{
    for (unsigned int i = 0; i < n_cargo.size(); i++)
        n_cargo[i]->kill_head(hd);
}

void cargo::motor_write(ostream& fout, int cargo_number)
{
    for (unsigned int i=0; i<n_cargo.size(); i++) {
        fout<<std::to_string(cargo_number)<<"\t"<<n_cargo[i]->write();
    }
}

/*void motor::update_angle()
{
    disp = rij_bc(BC, hx[1]-hx[0], hy[1]-hy[0], fov[0], fov[1], actin_network->get_delrx());
    mphi=atan2(disp[1],disp[0]);
}*/


/*array<double, 2> motor::boundary_check(int hd, double x, double y)
{
    return pos_bc(BC, actin_network->get_delrx(), dt, fov, {(x - hx[hd])/dt, (y - hy[hd])/dt}, {x, y});
}*/

/*string cargo::to_string()
{
    char buffer[1000];
    string out ="";

    sprintf(buffer, "\
            \nhead 0 position = (%f, %f)\t head 1 position=(%f,%f)\t angle = %f\
            \nstate = (%d, %d)\t f_index = (%d, %d)\t l_index = (%d, %d)\
            \nviscosity = %f\t max binding distance = %f\t stiffness = %f\t stall force = %f\t length = %f\
            \nkon = %f\t koff = %f\t kend = %f\t dt = %f\t temp = %f\t damp = %f\
            \nfov = (%f, %f)\t distance from end of link = (%f, %f)\
            shear = %f\t tension = (%f, %f)\n",
            hx[0], hy[0], hx[1], hy[1], mphi,
            state[0],  state[1], f_index[0],  f_index[1], l_index[0],  l_index[1],
            vs, max_bind_dist, mk, stall_force, mld,
            kon, koff, kend, dt, temperature, damp,
            fov[0],  fov[1], pos_a_end[0], pos_a_end[1], shear, force[0], force[1]);
    return buffer;
}*/


/*string cargo::write()
{
    return "\n" + std::to_string(hx[0]) + "\t" + std::to_string(hy[0])
        +  "\t" + std::to_string(disp[0]) + "\t" + std::to_string(disp[1])
        +  "\t" + std::to_string(f_index[0]) + "\t" + std::to_string(f_index[1])
        +  "\t" + std::to_string(l_index[0]) + "\t" + std::to_string(l_index[1]);
}*/

