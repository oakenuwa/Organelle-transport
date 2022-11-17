/*------------------------------------------------------------------
 motor_ensemble.cpp : container class for motors
 
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
#include "motor_ensemble.h"
#include "cargo.h"

#include "filament_ensemble.h"

//motor_ensemble class

motor_ensemble::motor_ensemble(double mdensity, int nCargo, array<double, 2> myfov, double delta_t, double temp,
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
    
    ke = 0;
    pe = 0;

    int total_mot = int(ceil(mdensity*fov[0]*fov[1]));
    if(nCargo>0){
        int nm = total_mot/nCargo;
        cout<<"\nDEBUG: Number of motors:"<<total_mot<<"\n";
        cout<<"\nDEBUG: Number of cargos:"<<nCargo<<"\n";

        //double alpha = 1, motorx, motory, mang;
        //array<double, 3> motor_pos;
        for (int i=0; i< nCargo; i++) {
            n_motors.push_back(new cargo(nm, fov, delta_t, temp, mlen, f_network,
                        v0, stiffness, max_ext_ratio, ron, roff, rend, fstall, rcut, vis, positions, BC));
            
        }
    }
}

motor_ensemble::motor_ensemble(double mdensity, array<double, 2> myfov, double delta_t, double temp,
      double mlen, filament_ensemble * network, double v0, double stiffness, double max_ext_ratio,
      double ron, double roff, double rend,
      double fstall, double rcut,
      double vis, vector<array<double,3> > positions, string BC) {

    fov = myfov;
    mld =mlen;
    gamma = 0;
    tMove=0;//10;
    f_network=network;
    v = v0;

    ke = 0;
    pe = 0;

    int nm = int(ceil(mdensity*fov[0]*fov[1]));
    cout<<"\nDEBUG: Number of motors:"<<nm<<"\n";

    double alpha = 1, motorx, motory, mang;
    array<double, 3> motor_pos;
    for (int i=0; i< nm; i++) {
      
      if ((unsigned int)i < positions.size()){
          motorx = positions[i][0];
          motory = positions[i][1];
          mang   = positions[i][2];
      }else{
          double theta = rng(0,2*pi);
          double d = rng(0,1);
          double len = sqrt(d)*11.0;
          motorx = len*cos(theta)-mld;//rng(-0.5*(fov[0]*alpha-mld),0.5*(fov[0]*alpha-mld));
          motory = len*sin(theta)-mld;//rng(-0.5*(fov[1]*alpha-mld),0.5*(fov[1]*alpha-mld));
          mang   = rng(0,2*pi);
      }
      motor_pos = {motorx, motory, mang};

      n_xlinks.push_back(new motor( motor_pos, mld, f_network,{0, 0}, {-1,-1}, {-1,-1}, fov, delta_t, temp,
                  v0, stiffness, max_ext_ratio, ron, roff, rend, fstall, rcut, vis, BC));
      
    }
}


motor_ensemble::motor_ensemble(vector<vector<double> > motors, array<double, 2> myfov, double delta_t, double temp, 
        double mlen, filament_ensemble * network, double v0, double stiffness, double max_ext_ratio, 
        double ron, double roff, double rend, 
        double fstall, double rcut,
        double vis, string BC) {
    
    fov = myfov;
    mld = mlen;
    gamma = 0;
    tMove = 0;
    f_network=network;
    v = v0;

    ke = 0;
    pe = 0;
    
    int nm = motors.size();
    cout<<"\nDEBUG: Number of motors:"<<nm<<"\n";

    array<double, 4> motor_pos;
    array<int, 2> f_index, l_index, state;

    for (int i=0; i< nm; i++) {
        
        motor_pos = {motors[i][0], motors[i][1], motors[i][2], motors[i][3]};
        
        f_index = {int(motors[i][4]), int(motors[i][5])};
        l_index = {int(motors[i][6]), int(motors[i][7])};

        state = {f_index[0] == -1 && l_index[0] == -1 ? 0 : 1, f_index[1] == -1 && l_index[1] == -1 ? 0 : 1};  

        n_xlinks.push_back(new motor( motor_pos, mld, f_network, state, f_index, l_index, fov, delta_t, temp,
                    v0, stiffness, max_ext_ratio, ron, roff, rend, fstall, rcut, vis, BC));
    }

    this->update_energies();
}

motor_ensemble::motor_ensemble(vector<vector<double> > motors, int nCargo, array<double, 2> myfov, double delta_t, double temp,
        double mlen, filament_ensemble * network, double v0, double stiffness, double max_ext_ratio,
        double ron, double roff, double rend,
        double fstall, double rcut,
        double vis, string BC) {
    
    fov = myfov;
    mld = mlen;
    gamma = 0;
    tMove = 0;
    f_network=network;
    v = v0;

    ke = 0;
    pe = 0;
    
    int nm = nCargo; // motors.size();
    cout<<"\nDEBUG: Number of cargos:"<<nm<<"\n";
    

    array<double, 4> motor_pos;
    array<int, 2> f_index, l_index, state;

    for (int i=0; i< nm; i++) {
        n_motors.push_back(new cargo(fov, delta_t, temp, f_network, vis, BC));
       // cout<<"did this work?"<<endl;
    }
    int totMot = motors.size();
    cout<<"\nDEBUG: Number of motors:"<<totMot<<"\n";
    for(int j=0; j<totMot; j++){
        int cargo_id = motors[j][0];
     //   cout<<"cargo_id: "<<cargo_id<<endl;
        motor_pos = {motors[j][1], motors[j][2], motors[j][3], motors[j][4]};
     //   cout<<cargo_id<<" "<<motors[j][1]<<" "<<motors[j][2]<<endl;
        f_index = {int(motors[j][5]), int(motors[j][6])};
        l_index = {int(motors[j][7]), int(motors[j][8])};

        state = {f_index[0] == -1 && l_index[0] == -1 ? 0 : 1, f_index[1] == -1 && l_index[1] == -1 ? 0 : 1};
        motor * q = new motor(motor_pos, mld, f_network, state, f_index, l_index, fov, delta_t, temp,
                      v0, stiffness, max_ext_ratio, ron, roff, rend, fstall, rcut, vis, BC, 1);
        n_motors[cargo_id]->add_motor(q);
    }
    
   // for(int j=0; j<nm; j++){
   //     for(int q=0; q<totMot/nm; q++)
   //     {
   //         array<double,2> p = n_motors[j]->get_hx(q);
   //         array<double,2> v = n_motors[j]->get_hy(q);
   //         cout<<"cargo: "<<j<<", motor: "<<q<<": "<<p[0]<<" "<<v[0]<<endl;
   //     }
   // }
    
    this->update_energies();
}


motor_ensemble::~motor_ensemble( ){ 
    cout<<"DELETING MOTOR ENSEMBLE\n";
    int s = n_motors.size();
    int q = n_xlinks.size();
    for (int i = 0; i < q; i++){
        delete n_xlinks[i];
    }
    
    for (int j = 0; j < s; j++){
        delete n_motors[j];
    }
    n_motors.clear();
    n_xlinks.clear();
};


int motor_ensemble::get_nmotors( ){ 
    return n_motors.size();
}


void motor_ensemble::kill_heads(int hd){
    for (unsigned int i = 0; i < n_xlinks.size(); i++)
        n_xlinks[i]->kill_head(hd);
}

//check if any motors attached to filaments that no longer exist; 
// if they do, detach them
// Worst case scenario, this is a O(m*n*p) function,
// where m is the number of filaments
//       n is the number of rods per filament
//       p is the number of motors
// However: we don't expect to fracture often, 
// so this loop should rarely if ever be accessed.
    


void motor_ensemble::check_broken_filaments()
{
    vector<int> broken_filaments = f_network->get_broken();
    array<int, 2> f_index;
    
    for (unsigned int i = 0; i < broken_filaments.size(); i++){
        
        for(unsigned int j = 0; j < n_xlinks.size(); j++){
            
            f_index = n_xlinks[j]->get_f_index();

            if(f_index[0] == broken_filaments[i]){
                n_xlinks[j]->detach_head_without_moving(0);
                //cout<<"\nDEBUG: detaching head 0 of motor "<<j<<" from broken filament "<<broken_filaments[i];
            }

            if(f_index[1] == broken_filaments[i]){
                n_xlinks[j]->detach_head_without_moving(1);
                //cout<<"\nDEBUG: detaching head 1 of motor "<<j<<" from broken filament "<<broken_filaments[i];
            }
        }
    }


}


void motor_ensemble::motor_walk(double t)
{

    this->check_broken_filaments();
    int nxlinks_sz = int(n_xlinks.size());
    bool attached;
    //#pragma omp parallel for
    
    for (int i=0; i<nxlinks_sz; i++) {
        
    //    if(i==0) cout<<"\nDEBUG: motor_walk: using "<<omp_get_num_threads()<<" cores";
        array<int, 2> s = n_xlinks[i]->get_states();
        
        if (t >= tMove){
            
            if (s[0] == 1)         
                n_xlinks[i]->step_onehead(0);
            else if (s[0] == 0)    
            {
                //cout<<"head zero about to attach"<<endl;
                attached = n_xlinks[i]->attach(0);
                if (!attached)
                    n_xlinks[i]->brownian_relax(0);
            }
            
            if (s[1] == 1)         
                n_xlinks[i]->step_onehead(1);
            else if (s[1] == 0)
            {
                //cout<<"head one about to attach"<<endl;
                attached = n_xlinks[i]->attach(1);
                if (!attached) 
                    n_xlinks[i]->brownian_relax(1);
            }

            n_xlinks[i]->update_angle();
            n_xlinks[i]->update_force();
            //n_motors[i]->update_force_fraenkel_fene();
            n_xlinks[i]->actin_update();
        }
    
    }
    
    this->update_energies();
    
}

void motor_ensemble::cargo_walk(double t)
{
    int nmotors_sz = int(n_motors.size());
    //this->check_broken_filaments();
    for (int j = 0; j<nmotors_sz; j++)
    {
        cout<<"Cargo number: "<<j<<endl;
        n_motors[j]->motor_walk(t);
    }

}
/* Used for static, contstantly attached, motors -- ASSUMES both heads are ALWAYS attached */

void motor_ensemble::motor_update()
{

    this->check_broken_filaments();
    int nxlinks_sz = int(n_xlinks.size());
    //#pragma omp parallel for
    
    for (int i=0; i<nxlinks_sz; i++) {
       
            n_xlinks[i]->update_position_attached(0);
            n_xlinks[i]->update_position_attached(1);
            n_xlinks[i]->update_angle();
            n_xlinks[i]->update_force();
            //n_motors[i]->update_force_fraenkel_fene();
            n_xlinks[i]->actin_update();
    
    }
    this->update_energies();
    
}

void motor_ensemble::xlink_write(ostream& fout)
{
    for (unsigned int i=0; i<n_xlinks.size(); i++) {
        fout<<n_xlinks[i]->write();
    }
}

void motor_ensemble::motor_write(ostream& fout)
{
    for (unsigned int i=0; i<n_motors.size(); i++) {
        n_motors[i]->motor_write(fout,i);
    }
}

void motor_ensemble::add_motor(motor * m)
{
    n_xlinks.push_back(m);
}


void motor_ensemble::set_shear(double g)
{
    for (unsigned int i=0; i<n_xlinks.size(); i++)
        n_xlinks[i]->set_shear(g);
    
    gamma = g;
}

 
void motor_ensemble::update_energies()
{
    ke = 0;
    pe = 0;
    for (unsigned int m = 0; m < n_xlinks.size(); m++)
    {
        ke += n_xlinks[m]->get_kinetic_energy();
        pe += n_xlinks[m]->get_stretching_energy();
        //pe += n_motors[m]->get_stretching_energy_fene();
    }
}

 
double motor_ensemble::get_potential_energy(){
    return pe;
}

 
void motor_ensemble::print_ensemble_thermo(){
    cout<<"\nAll Motors\t:\tKE = "<<ke<<"\tPEs = "<<pe<<"\tPEb = "<<0<<"\tTE = "<<(ke+pe);
}
