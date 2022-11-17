/*
 * motor.h
 *  
 *
 *  Created by Shiladitya Banerjee on 9/3/13.
 *  Copyright 2013 University of Chicago. All rights reserved.
 *
 */

//=====================================
//include guard
#ifndef __MOTOR_H_INCLUDED__
#define __MOTOR_H_INCLUDED__

//=====================================
// forward declared dependencies
class filament_ensemble;

//=====================================
//included dependences
#include "globals.h"

//motor class
class motor 
{
    public:

        motor(array<double, 3> pos, double mlen, filament_ensemble* network, 
                array<int, 2> mystate, array<int, 2> myfindex, array<int, 2> myrindex,
                array<double, 2> myfov, double delta_t, double v0, double temp, double stiffness, double max_ext_ratio, 
                double ron, double roff, double rend, 
                double fstall, double rcut,
                double vis, string BC);
    
    motor(array<double, 3> pos, double mlen, filament_ensemble* network,
            array<int, 2> mystate, array<int, 2> myfindex, array<int, 2> myrindex,
            array<double, 2> myfov, double delta_t, double v0, double temp, double stiffness, double max_ext_ratio,
            double ron, double roff, double rend,
            double fstall, double rcut,
            double vis, string BC, int indicator);
        
        motor(array<double, 4> pos, double mlen, filament_ensemble* network, 
                array<int, 2> mystate, array<int, 2> myfindex, array<int, 2> myrindex,
                array<double, 2> myfov, double delta_t, double v0, double temp, double stiffness, double max_ext_ratio,
                double ron, double roff, double rend, 
                double fstall, double rcut,
                double vis, string BC);
    
    motor(array<double, 4> pos, double mlen, filament_ensemble* network,
            array<int, 2> mystate, array<int, 2> myfindex, array<int, 2> myrindex,
            array<double, 2> myfov, double delta_t, double v0, double temp, double stiffness, double max_ext_ratio,
            double ron, double roff, double rend,
            double fstall, double rcut,
            double vis, string BC, int indicator);
        
        motor();
        
        ~motor();

        string get_BC();

        bool allowed_bind( int hd, array<int, 2> fl_idx);
        
        bool attach( int hd);
    
        bool attach_motors( int hd);

        void relax_head( int hd);

        void kill_head( int hd);
        
        void update_position_attached(int hd);
        
        void update_force();
        
        void update_force_fraenkel_fene();
        
        void update_angle();
        
        void brownian_relax(int hd);

        array<double, 2> boundary_check(int i,  double vx, double vy);

        void step_onehead( int hd);
    
        void step_onehead_motors( int hd);

        void step_twoheads();

        void actin_update_hd(int hd, array<double, 2> f);
        
        void actin_update();
    
        void actin_update_angular();

        void update_shape();
        
        void set_shear(double g);
    
        void set_hx(int hd, double posx);
        
        void set_hy(int hd, double posy);
    
        array<int,2> get_f_index();
        
        array<int,2> get_l_index();
        
        array<double,2> get_pos_a_end();
        
        array<double,2> get_force();
        
        array<int, 2> get_states();
        
        array<double, 2> get_hx();

        array<double, 2> get_hy();

        void detach_head(int hd, array<double, 2> pos);

        void detach_head_without_moving(int hd);
        
        void update_pos_a_end(int hd, double pos);

        inline void reflect(double t, double gamma, double x1, double x2, double y1, double y2);
        
        inline void periodic(double t, double gamma, double x1, double x2, double y1, double y2);

        double get_stretching_energy();
        
        double get_stretching_energy_fene();

        double get_kinetic_energy();
        
        double metropolis_prob(int hd, array<int, 2> fl_idx, array<double, 2> newpos, double maxprob);
    
       double metropolis_prob_angle(int hd, array<int, 2> fl_idx, array<double, 2> newpos, double maxprob, double angle);

        array<double, 2> generate_off_pos(int hd);

        string to_string();
        
        string write();
    
    public:

        double mphi,mld, vs, stall_force, max_bind_dist, mk, kon, koff, kend, dt, temperature, 
               damp, shear, max_ext, eps_ext, kinetic_energy, bd_prefactor, tension;
        
        array<double,2> hx, hy, pos_a_end, fov, prv_rnd_x, prv_rnd_y, force, disp;

        array<array<double, 2>, 2> ldir_bind, bind_disp;

        array<int,2> state, f_index, l_index;
        
        map<vector<int>, double> dist;
        array<bool, 2> at_barbed_end;
        
        string BC;
        
        filament_ensemble* actin_network;
    
        int indicator;
        
};

#endif
