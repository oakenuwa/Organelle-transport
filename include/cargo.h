//
//  cargo.hpp
//  
//
//  Created by Oghosa Akenuwa on 9/27/21.
//

#ifndef __CARGO_H_INCLUDED__
#define __CARGO_H_INCLUDED__

//=====================================
// forward declared dependencies
//class filament_ensemble;
class motor;
//=====================================
//included dependences
#include "motor.h"
//#include "globals.h"
//#include "motor.h";
#include "vector"

//motor ensemble class

class cargo
{
    public:
    
        cargo(int mpercargo, array<double, 2> myfov, double delta_t, double temp,
            double mlen, filament_ensemble * network, double v0, double stiffness, double max_ext_ratio,
            double ron, double roff, double rend,
            double fstall, double rcut,
            double vis, vector<array<double,3> > positions, string BC);
    
        cargo(array<double, 2> myfov, double delta_t, double temp,
              filament_ensemble * network,
              double vis, string BC);
    
        ~cargo();
    
//    void brownian_relax();
    
    array<int, 2> get_states(int);
    
    array<double, 2> get_hx(int);

    array<double, 2> get_hy(int);
    
    string get_BC(int);
    

//        int get_nmotors();

        
        void motor_walk(double t);

//        void motor_update();
        
//        void update_energies();
        
//        double get_potential_energy();

        void motor_write(ostream& fout, int cargo_number);

 //       void print_ensemble_thermo();
        
 //       void motor_tension(ofstream& fout);

        void add_motor(motor * m);

//        void set_shear(double g);
//        vector<motor *> n_cargo_test;
    
    private:

        double mld, gamma, tMove, cmotorx, cmotory, cdamp, cbd_prefactor, cdt, ctemperature;
        double ke, pe, v, c_prv_rnd_x, c_prv_rnd_y;
        string bc;
        void kill_heads(int i);
        void check_broken_filaments();
        array<double, 2> cargo_boundary_check(double x, double y);
        array<double, 2> fov;
        filament_ensemble *f_network;
        vector<motor *> n_cargo;
};

#endif
