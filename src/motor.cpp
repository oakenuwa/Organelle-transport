/*------------------------------------------------------------------
 motor.cpp : object describing a motor or crosslinker
 
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
#include "filament_ensemble.h"
//#include "actin.h"

//motor class

motor::motor(){}

motor::motor( array<double, 3> pos, 
        double mlen, filament_ensemble * network, 
        array<int, 2> mystate, 
        array<int, 2> myfindex, 
        array<int, 2> mylindex,
        array<double, 2> myfov,
        double delta_t, 
        double temp,
        double v0, 
        double stiffness, 
        double max_ext_ratio, 
        double ron, double roff, double rend, 
        double fstall, double rcut,
        double vis, string bc) {
    
    vs          = v0;
    mk          = stiffness;//rng(10,100); 
    
    stall_force   = fstall;
    temperature   = temp;
    
    max_bind_dist = rcut;

    mld         = mlen;
    dt          = delta_t;
    kon         = ron*dt;
    koff        = roff*dt;
    kend        = rend*dt;
    mphi        = pos[2];
    state       = mystate;
    f_index     = myfindex; //filament index for each head
    l_index     = mylindex; //link index for each head
    
    fov         = myfov;
    BC          = bc; 
    actin_network = network;
    damp        = (6*pi*vis*mld);
    bd_prefactor= sqrt(temperature/(2*damp*dt)); 

    /****for FENE motors******/
    max_ext     = max_ext_ratio*mlen;
    eps_ext     = 0.01*max_ext;
    /*************************/

    shear       = 0;
    tension     = 0;
    force       = {0,0}; // force on the spring  
    kinetic_energy = 0; //assume m = 1
    
    array<double, 2> posH0 = boundary_check(0, pos[0]-0.5*mld*cos(mphi), pos[1]-0.5*mld*sin(mphi));
    array<double, 2> posH1 = boundary_check(1, pos[0]+0.5*mld*cos(mphi), pos[1]+0.5*mld*sin(mphi));
    hx[0] = posH0[0];
    hy[0] = posH0[1];
    hx[1] = posH1[0];
    hy[1] = posH1[1];
    
    disp = rij_bc(BC, hx[1]-hx[0], hy[1]-hy[0], fov[0], fov[1], actin_network->get_delrx()); 
    
    pos_a_end = {0, 0}; // pos_a_end = distance from pointy end -- by default 0
                        // i.e., if l_index[hd] = j, then pos_a_end[hd] is the distance to the "j+1"th actin
    
    ldir_bind[0] = {0,0};
    ldir_bind[1] = {0,0};

    bind_disp[0] = {0,0};
    bind_disp[1]=  {0,0};

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

motor::motor( array<double, 3> pos,
        double mlen, filament_ensemble * network,
        array<int, 2> mystate,
        array<int, 2> myfindex,
        array<int, 2> mylindex,
        array<double, 2> myfov,
        double delta_t,
        double temp,
        double v0,
        double stiffness,
        double max_ext_ratio,
        double ron, double roff, double rend,
        double fstall, double rcut,
        double vis, string bc, int indicator) {
    
    vs          = v0;
    mk          = stiffness;//rng(10,100);
    
    stall_force   = fstall;
    temperature   = temp;
    
    max_bind_dist = rcut;

    mld         = mlen;
    dt          = delta_t;
    kon         = ron*dt;
    koff        = roff*dt;
    kend        = rend*dt;
    mphi        = pos[2];
    state       = mystate;
    f_index     = myfindex; //filament index for each head
    l_index     = mylindex; //link index for each head
    
    fov         = myfov;
    BC          = bc;
    actin_network = network;
    damp        = (6*pi*vis*mld);
    bd_prefactor= sqrt(temperature/(2*damp*dt));

    /****for FENE motors******/
    max_ext     = max_ext_ratio*mlen;
    eps_ext     = 0.01*max_ext;
    /*************************/

    shear       = 0;
    tension     = 0;
    force       = {0,0}; // force on the spring
    kinetic_energy = 0; //assume m = 1
    
    array<double, 2> posH0 = boundary_check(0, pos[0], pos[1]);
    array<double, 2> posH1 = boundary_check(1, pos[0]+mld*cos(mphi), pos[1]+mld*sin(mphi));
    
    hx[0] = posH0[0];
    hy[0] = posH0[1];
    hx[1] = posH1[0];
    hy[1] = posH1[1];
    
    disp = rij_bc(BC, hx[1]-hx[0], hy[1]-hy[0], fov[0], fov[1], actin_network->get_delrx());
    
    pos_a_end = {0, 0}; // pos_a_end = distance from pointy end -- by default 0
                        // i.e., if l_index[hd] = j, then pos_a_end[hd] is the distance to the "j+1"th actin
    
    ldir_bind[0] = {0,0};
    ldir_bind[1] = {0,0};

    bind_disp[0] = {0,0};
    bind_disp[1]=  {0,0};

    at_barbed_end = {false, false};
    
    cout<<"posH1: "<<mld*cos(mphi)<<", posH1: "<<mld*sin(mphi)<<" posH1: "<<posH1[0]<<", posH1: "<<posH1[1]<<endl;
    
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

motor::motor( array<double, 4> pos, 
        double mlen, filament_ensemble * network, 
        array<int, 2> mystate, 
        array<int, 2> myfindex, 
        array<int, 2> mylindex,
        array<double, 2> myfov,
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
    
    /********for FENE springs*********/ 
    max_ext     = max_ext_ratio*mlen;
    eps_ext     = 0.01*max_ext;
    /********************************/

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

motor::motor( array<double, 4> pos,
        double mlen, filament_ensemble * network,
        array<int, 2> mystate,
        array<int, 2> myfindex,
        array<int, 2> mylindex,
        array<double, 2> myfov,
        double delta_t,
        double temp,
        double v0,
        double stiffness,
        double max_ext_ratio,
        double ron, double roff, double rend,
        double fstall, double rcut,
        double vis, string bc, int indicator) {
    
    vs          = v0;
    mk          = stiffness;//rng(10,100);
    
    stall_force   = fstall;
    temperature   = temp;
    
    max_bind_dist = rcut;

    mld         = mlen;
    dt          = delta_t;
    kon         = ron*dt;
    koff        = roff*dt;
    kend        = rend*dt;
    mphi        = pos[2];
    state       = mystate;
    f_index     = myfindex; //filament index for each head
    l_index     = mylindex; //link index for each head
    
    fov         = myfov;
    BC          = bc;
    actin_network = network;
    damp        = (6*pi*vis*mld);
    bd_prefactor= sqrt(temperature/(2*damp*dt));

    /****for FENE motors******/
    max_ext     = max_ext_ratio*mlen;
    eps_ext     = 0.01*max_ext;
    /*************************/

    shear       = 0;
    tension     = 0;
    force       = {0,0}; // force on the spring
    kinetic_energy = 0; //assume m = 1
    
    array<double, 2> posH0 = boundary_check(0, pos[0], pos[1]);
    array<double, 2> posH1 = boundary_check(1, pos[0]+pos[2], pos[1]+pos[3]);
    hx[0] = posH0[0];
    hy[0] = posH0[1];
    hx[1] = posH1[0];
    hy[1] = posH1[1];
    
    disp = rij_bc(BC, hx[1]-hx[0], hy[1]-hy[0], fov[0], fov[1], actin_network->get_delrx());
    
    pos_a_end = {0, 0}; // pos_a_end = distance from pointy end -- by default 0
                        // i.e., if l_index[hd] = j, then pos_a_end[hd] is the distance to the "j+1"th actin
    
    ldir_bind[0] = {0,0};
    ldir_bind[1] = {0,0};

    bind_disp[0] = {0,0};
    bind_disp[1]=  {0,0};

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


 motor::~motor(){};

//return motor state with a given head number

array<int, 2> motor::get_states() 
{
    return state;
}

array<double, 2> motor::get_hx()
{
    return hx;
}


array<double, 2> motor::get_hy()
{
    return hy;
}


string motor::get_BC()
{
    return BC;
}


void motor::set_shear(double gamma)
{
    shear = gamma;
}

void motor::set_hx(int hd, double posx)
{
    hx[hd] = posx;
}

void motor::set_hy(int hd, double posy)
{
    hy[hd] = posy;
}

//metropolis algorithm with rate constant
//NOTE: while fl_idx doesn't matter for this xlink implementation, it does for "spacers"
double motor::metropolis_prob(int hd, array<int, 2> fl_idx, array<double, 2> newpos, double maxprob)
{
    //array<double,2> force = this->get_force();
    double prob = maxprob; //*exp(sqrt(force[0]*force[0] + force[1]*force[1])/5);
    double stretch  = dist_bc(BC, newpos[0] - hx[pr(hd)], newpos[1] - hy[pr(hd)], fov[0], fov[1], actin_network->get_delrx()) - mld; 
    double delE = 0.5*mk*stretch*stretch - this->get_stretching_energy();

    if( delE > 0 )
        prob *= exp(-delE/temperature);
    
    return prob;
}

double motor::metropolis_prob_angle(int hd, array<int, 2> fl_idx, array<double, 2> newpos, double maxprob, double del_angle)
{

    double prob = maxprob;
    double stretch  = dist_bc(BC, newpos[0] - hx[pr(hd)], newpos[1] - hy[pr(hd)], fov[0], fov[1], actin_network->get_delrx()) - mld;
    double delE = 0.5*mk*stretch*stretch - this->get_stretching_energy();
    double del_E = delE + del_angle;
    if( del_E > 0)
    {
        prob *= exp(-del_E/temperature);
    }
    
    return prob;
}

bool motor::allowed_bind(int hd, array<int, 2> fl_idx){
    return (f_index[pr(hd)] != fl_idx[0] || l_index[pr(hd)] != fl_idx[1]);
}

//check for attachment of unbound heads given head index (0 for head 1, and 1 for head 2)
bool motor::attach(int hd)
{
    
    double not_off_prob = 0;
    double mf_rand = rng(0,1.0);
    array<double, 2> intPoint;

//    set<pair<double, array<int, 2> > > dist_sorted = actin_network->get_dist_all(hx[hd], hy[hd]);//if not using neighbor lists
    set<pair<double, array<int, 2> > > dist_sorted = actin_network->get_dist(hx[hd], hy[hd]);
    //if(dist_sorted.empty())
    //{//cout<<"No links found in sight"<<endl;}
    if(!dist_sorted.empty()){

        for (set<pair<double, array<int, 2>>>::iterator it=dist_sorted.begin(); it!=dist_sorted.end(); ++it)
        {
            //cout<<"f_index:l_index: "<<endl; //(it->second).at(0)<<" "<<(it->second).at(1)<<endl;
            if (it->first > max_bind_dist) //since it's sorted, all the others will be farther than max_bind_dist too
                break;

            //head can't bind to the same filament link the other head is bound to
        //else if(!(f_index[pr(hd)]==(it->second).at(0) && l_index[pr(hd)]==(it->second).at(1))) {
            
            
            else if(allowed_bind(hd, it->second))
            {
                //cout<<"Test for motor attach: "<<hd<<endl;
                intPoint = actin_network->get_filament((it->second).at(0))->get_link((it->second).at(1))->get_intpoint();
                not_off_prob += metropolis_prob(hd, it->second, intPoint, kon);

                if (mf_rand < not_off_prob && state[pr(hd)] == 0)
                {
                    //cout<< "single xlinker head attachment came here"<<endl;
                    //if (state[pr(hd)] == 0)
                    //{
                        //cout<<"this worked too"<<endl;
                        state[hd] = 1;
                        f_index[hd] = (it->second).at(0);
                        l_index[hd] = (it->second).at(1);

                        //record displacement of head and orientation of link for future unbinding move
                        ldir_bind[hd] = actin_network->get_direction(f_index[hd], l_index[hd]);
                        bind_disp[hd] = rij_bc(BC, intPoint[0]-hx[hd], intPoint[1]-hy[hd], fov[0], fov[1], actin_network->get_delrx());

                        //update head position
                        hx[hd] = intPoint[0];
                        hy[hd] = intPoint[1];

                        //update relative head position
                        pos_a_end[hd]=dist_bc(BC, actin_network->get_end(f_index[hd], l_index[hd])[0] - hx[hd],
                                              actin_network->get_end(f_index[hd], l_index[hd])[1] - hy[hd], fov[0], fov[1],
                                              actin_network->get_delrx());

                        //(even if its at the barbed end upon binding, could have negative velocity, so always set this to false, until it steps)
                        at_barbed_end[hd] = false;

                        return true;
                    }

                    else if (state[pr(hd)] == 1)
                    {
                        double new_not_off_prob = 0.0;
                        //cout<<"new_metropolis_factor for xlinker head attachment instance"<<endl;
                        int x = f_index[pr(hd)];
                        int y = l_index[pr(hd)];
                        int z = (it->second).at(0);
                        int w = (it->second).at(1);
                        array<double, 2> disp1 = actin_network->get_filament(x)->get_link(y)->get_neg_disp();
                        array<double, 2> disp2 = actin_network->get_filament(z)->get_link(w)->get_neg_disp();
                        //double dotProduct = dot(disp1, disp2);
                        
                        
                        double rsq1  = disp1[0]*disp1[0] + disp1[1]*disp1[1];
                        double r1    = sqrt(rsq1);
                        
                        // 2nd bond
                        //delr2 = links[n+1]->get_disp();
                        double rsq2  = disp2[0]*disp2[0] + disp2[1]*disp2[1];
                        double r2    = sqrt(rsq2);

                        // angle (cos and sin)
                        double c = (disp1[0]*disp2[0] + disp1[1]*disp2[1]) / (r1*r2);
                       
                        double a   = 0;//0.011 * (acos(c) - 0) * (acos(c) - 0)/2;
                        new_not_off_prob += metropolis_prob_angle(hd, it->second, intPoint, kon, a);

                        //if (dotProduct < 0.9)
                        //{
                            //cout<<"this worked as well"<<endl;
                        //    break;
                        //}

                        if (mf_rand < new_not_off_prob)                                            //(dotProduct >= 0.9 && dotProduct <= 1.0)
                        {
                            //cout<<"double crosslinking instance"<<endl;
                            //cout<<"\nangle between links: "<< acos(c)*180/pi<<endl;
                            
                            //cout<<"\n"<< acos(c)*180/pi;
                            
                            //update state
                            state[hd] = 1;
                            f_index[hd] = (it->second).at(0);
                            l_index[hd] = (it->second).at(1);

                            //record displacement of head and orientation of link for future unbinding move
                            ldir_bind[hd] = actin_network->get_direction(f_index[hd], l_index[hd]);
                            bind_disp[hd] = rij_bc(BC, intPoint[0]-hx[hd], intPoint[1]-hy[hd], fov[0], fov[1], actin_network->get_delrx());

                            //update head position
                            hx[hd] = intPoint[0];
                            hy[hd] = intPoint[1];

                            //update relative head position
                            pos_a_end[hd]=dist_bc(BC, actin_network->get_end(f_index[hd], l_index[hd])[0] - hx[hd],
                                                  actin_network->get_end(f_index[hd], l_index[hd])[1] - hy[hd], fov[0], fov[1],
                                                  actin_network->get_delrx());

                            //(even if its at the barbed end upon binding, could have negative velocity, so always set this to false, until it steps)
                            at_barbed_end[hd] = false;

                            return true;
                        }
                }
            }
        }
    }
    return false;
} 

bool motor::attach_motors(int hd)
{
    double not_off_prob = 0;
    double mf_rand = rng(0,1.0);
    array<double, 2> intPoint;
    
//    set<pair<double, array<int, 2> > > dist_sorted = actin_network->get_dist_all(hx[hd], hy[hd]);//if not using neighbor lists
    set<pair<double, array<int, 2> > > dist_sorted = actin_network->get_dist(hx[hd], hy[hd]);

    if(!dist_sorted.empty()){
        
        for (set<pair<double, array<int, 2>>>::iterator it=dist_sorted.begin(); it!=dist_sorted.end(); ++it)
        {
            if (it->first > max_bind_dist) //since it's sorted, all the others will be farther than max_bind_dist too
                break;

            //head can't bind to the same filament link the other head is bound to
 //           else if(!(f_index[pr(hd)]==(it->second).at(0) && l_index[pr(hd)]==(it->second).at(1))) {
            else if(allowed_bind(hd, it->second))
            {
                
                intPoint = actin_network->get_filament((it->second).at(0))->get_link((it->second).at(1))->get_intpoint();
                not_off_prob += metropolis_prob(hd, it->second, intPoint, kon);
                 
                if (mf_rand < not_off_prob)
                {
                    //cout<< "this came here"<<endl;
                        //cout<<"this worked too"<<endl;
                        state[hd] = 1;
                        f_index[hd] = (it->second).at(0);
                        l_index[hd] = (it->second).at(1);
                        
                        //record displacement of head and orientation of link for future unbinding move
                        ldir_bind[hd] = actin_network->get_direction(f_index[hd], l_index[hd]);
                        bind_disp[hd] = rij_bc(BC, intPoint[0]-hx[hd], intPoint[1]-hy[hd], fov[0], fov[1], actin_network->get_delrx());
                        
                        //update head position
                        hx[hd] = intPoint[0];
                        hy[hd] = intPoint[1];
                        
                        //update relative head position
                        pos_a_end[hd]=dist_bc(BC, actin_network->get_end(f_index[hd], l_index[hd])[0] - hx[hd],
                                              actin_network->get_end(f_index[hd], l_index[hd])[1] - hy[hd], fov[0], fov[1],
                                              actin_network->get_delrx());
                        
                        //(even if its at the barbed end upon binding, could have negative velocity, so always set this to false, until it steps)
                        at_barbed_end[hd] = false;
                        
                        return true;
                }
            }
        }
    }
    return false;
}

void motor::update_force()
{ 
    //force = {mk*(disp[0]-mld*cos(mphi)), mk*(disp[1]-mld*sin(mphi))};
    tension = mk*(hypot(disp[0], disp[1]) - mld);
    force = {tension*cos(mphi), tension*sin(mphi)};
}

/* Taken from hsieh, jain, larson, jcp 2006; eqn (5)
 * Adapted by placing a cutoff, similar to how it's done in LAMMPS src/bond_fene.cpp*/

void motor::update_force_fraenkel_fene()
{
    double ext = abs(mld - hypot(disp[0], disp[1]));
    double scaled_ext, mkp;
    
    if (max_ext - ext > eps_ext )
        scaled_ext = ext/max_ext;
    else
        scaled_ext = (max_ext - eps_ext)/max_ext;
    
    mkp = mk/(1-scaled_ext*scaled_ext);
    force = {mkp*(disp[0]-mld*cos(mphi)), mkp*(disp[1]-mld*sin(mphi))};

}


void motor::brownian_relax(int hd)
{
    array<double,2>forces = this->get_force();
    //double k = 1; //pN/um
    //double r = 11.2838; //um
    //double a = r; //semi-major axis
    //double b = r; //semi-minor axis
    //double posx = hx[hd];
    //double posy = hy[hd];
    //double m = posy/posx;
    //double lsqr = pow(posx,2)+pow(posy,2);
    //double Fx = 0;
    //double Fy = 0;
    //double ellip_check = pow(posx/a,2)+pow(posy/b,2);
    /*if (lsqr>=(r*r))
    {
        double l = sqrt(lsqr);
        double dl = l-r;
        Fx = -k*dl*(posx/l);
        Fy = -k*dl*(posy/l);
        //cout<<"confine"<<endl;

    
    
       /* if (ellip_check>=1)
        {
            double l = sqrt(lsqr);
            double l_act_ellip = sqrt((pow(a,2)*pow(b,2)*(1+pow(m,2))/(pow(b,2)+pow(a,2)*pow(m,2))));
            double dl = l-l_act_ellip;
            double Fx = -k*dl*(posx/l);
            double Fy = -k*dl*(posy/l);
            force[0]+=Fx;
            force[1]+=Fy;
            //cout<<"confine"<<endl;
       }
        double new_rnd_x= rng_n(0,1), new_rnd_y = rng_n(0,1);
        
        double vx =  (pow(-1,hd)*forces[0] + Fx) / damp + bd_prefactor*(new_rnd_x + prv_rnd_x[hd]);
        double vy =  (pow(-1,hd)*forces[1] + Fy)/ damp + bd_prefactor*(new_rnd_y + prv_rnd_y[hd]);
        kinetic_energy = vx*vx + vy*vy;
        array<double, 2> pos = boundary_check(hd, hx[hd] + vx*dt, hy[hd] + vy*dt);
        hx[hd] = pos[0];
        hy[hd] = pos[1];

        prv_rnd_x[hd] = new_rnd_x;
        prv_rnd_y[hd] = new_rnd_y;
    }*/
    
    //else
    //{
        double new_rnd_x= rng_n(0,1), new_rnd_y = rng_n(0,1);
        
        double vx =  (pow(-1,hd)*forces[0]) / damp + bd_prefactor*(new_rnd_x + prv_rnd_x[hd]);
        double vy =  (pow(-1,hd)*forces[1]) / damp + bd_prefactor*(new_rnd_y + prv_rnd_y[hd]);
        kinetic_energy = vx*vx + vy*vy;
        array<double, 2> pos = boundary_check(hd, hx[hd] + vx*dt, hy[hd] + vy*dt);
        hx[hd] = pos[0];
        hy[hd] = pos[1];

        prv_rnd_x[hd] = new_rnd_x;
        prv_rnd_y[hd] = new_rnd_y;
    //}

}


void motor::kill_head(int hd)
{
    state[hd] = -1;
}


void motor::relax_head(int hd)
{
    array<double, 2> newpos = boundary_check(hd, hx[pr(hd)] - pow(-1, hd)*mld*cos(mphi), hy[pr(hd)] - pow(-1, hd)*mld*sin(mphi));
    hx[hd] = newpos[0];
    hy[hd] = newpos[1];
}


void motor::update_angle()
{
    disp = rij_bc(BC, hx[1]-hx[0], hy[1]-hy[0], fov[0], fov[1], actin_network->get_delrx()); 
    mphi=atan2(disp[1],disp[0]);
}


array<double, 2> motor::boundary_check(int hd, double x, double y)
{
    return pos_bc(BC, actin_network->get_delrx(), dt, fov, {(x - hx[hd])/dt, (y - hy[hd])/dt}, {x, y});
}

array<double, 2> motor::generate_off_pos(int hd){
    
    array<double, 2> ldir = actin_network->get_direction(f_index[hd], l_index[hd]);
    double c = dot(  ldir, ldir_bind[hd]);
    double s = cross(ldir, ldir_bind[hd]);

    array<double, 2> bind_disp_rot = {bind_disp[hd][0]*c - bind_disp[hd][1]*s, bind_disp[hd][0]*s + bind_disp[hd][1]*c};

    return pos_bc(BC, actin_network->get_delrx(), dt, fov, 
            {-bind_disp_rot[0]/dt, -bind_disp_rot[1]/dt}, 
            {hx[hd] - bind_disp_rot[0], hy[hd] - bind_disp_rot[1]}
            ); 
    //array<double, 2> newpos = {hx[hd]-bind_disp_rot[0], hy[hd]-bind_disp_rot[1]};
    //return boundary_check(hd, newpos[0], newpos[1]);
} 


//stepping and detachment kinetics of a single bound head 
void motor::step_onehead(int hd)
{
    array<double, 2> hpos_new = generate_off_pos(hd);
    
    if (state[pr(hd)] == 0)
    {
        double a_ngle = 0;
        double off_prob = metropolis_prob_angle(hd, {0,0}, hpos_new, at_barbed_end[hd] ? kend : koff,a_ngle);
        if ( event(off_prob) ) this->detach_head(hd, hpos_new);
        else{

            //calculate motor velocity
            if (vs != 0 && !(at_barbed_end[hd])){
                double vm = vs;
                if (state[pr(hd)] != 0){
                    vm = my_velocity(vs,
                            pow(-1, hd)*dot(force, actin_network->get_direction(f_index[hd], l_index[hd])),
                            stall_force);
                }
                this->update_pos_a_end(hd, pos_a_end[hd]+dt*vm); // update relative position
            }
            if (state[hd] == 1) this->update_position_attached(hd);  // update absolute position

        }
    
    }
    
    else if (state[pr(hd)] != 0)
    {
//        //double off_prob = 0.0;
        //cout<<"\n for xlinker head detachment instance"<<endl;
        int x = f_index[pr(hd)];
        int y = l_index[pr(hd)];
        int z = f_index[hd];
        int w = l_index[hd];
        //cout<<"\nDEBUG: f_index0, l_index0, f_index1, l_index1 = "<<f_index[hd]<<" , "<<l_index[hd]<<" , "<<f_index[pr(hd)]<<" , "<<l_index[pr(hd)];
        
        array<double, 2> disp_1 = actin_network->get_filament(x)->get_link(y)->get_neg_disp();
        array<double, 2> disp_2 = actin_network->get_filament(z)->get_link(w)->get_neg_disp();
//        //double dotProduct = dot(disp1, disp2);
//
//
        double rsq_1  = disp_1[0]*disp_1[0] + disp_1[1]*disp_1[1];
        double r_1    = sqrt(rsq_1);
//
//        // 2nd bond
//        //delr2 = links[n+1]->get_disp();
        double rsq_2  = disp_2[0]*disp_2[0] + disp_2[1]*disp_2[1];
        double r_2    = sqrt(rsq_2);
//
//        // angle (cos and sin)
        double c = (disp_1[0]*disp_2[0] + disp_1[1]*disp_2[1]) / (r_1*r_2);
//
        double a   = -0.011 * (acos(c) - 0) * (acos(c) - 0)/2;
//        //double off_prob = metropolis_prob_angle(hd, it->second, intPoint, kon, a);
//
//        //double off_prob = metropolis_prob_angle(hd, {0,0}, hpos_new, at_barbed_end[hd] ? kend : koff,-a);
        double new_off_prob = metropolis_prob_angle(hd, {0,0}, hpos_new, at_barbed_end[hd] ? kend : koff,a);
//        //cout<<"off_prob: "<<off_prob<<endl;
        //cout<<"\nDEBUG: at barbed end? : "<<at_barbed_end[hd]<<"; off_prob = "<<new_off_prob;
//
//        // attempt detachment
        if ( event(new_off_prob) ) this->detach_head(hd, hpos_new);
//
        else{

            //calculate motor velocity
            if (vs != 0 && !(at_barbed_end[hd])){
                double vm = vs;
                if (state[pr(hd)] != 0)
                {
                    vm = my_velocity(vs,
                            pow(-1, hd)*dot(force, actin_network->get_direction(f_index[hd], l_index[hd])),
                            stall_force);
                }
                this->update_pos_a_end(hd, pos_a_end[hd]+dt*vm); // update relative position
            }
            if (state[hd] == 1) this->update_position_attached(hd);  // update absolute position

        }
    }

}

void motor::step_onehead_motors(int hd)
{

    array<double, 2> hpos_new = generate_off_pos(hd);
    double off_prob = metropolis_prob(hd, {0,0}, hpos_new, at_barbed_end[hd] ? kend : koff);
    
    //cout<<"\nDEBUG: at barbed end? : "<<at_barbed_end[hd]<<"; off_prob = "<<off_prob;
    // attempt detachment
    if ( event(off_prob) ) this->detach_head(hd, hpos_new);
    else{

        //calculate motor velocity
        if (vs != 0 && !(at_barbed_end[hd])){
            double vm = vs;
            if (state[pr(hd)] != 0){
                vm = my_velocity(vs,
                        pow(-1, hd)*dot(force, actin_network->get_direction(f_index[hd], l_index[hd])),
                        stall_force);
            }
            this->update_pos_a_end(hd, pos_a_end[hd]+dt*vm); // update relative position
        }
        if (state[hd] == 1) this->update_position_attached(hd);  // update absolute position
        
    }
}


void motor::update_pos_a_end(int hd, double pos)
{
//    cout<<"\nDEBUG: new pos = "<<pos;
    double link_length = actin_network->get_llength(f_index[hd],l_index[hd]);
    if (pos >= link_length) { // "passed" the link
        if (l_index[hd] == 0){ // the barbed end of the filament
            at_barbed_end[hd] = true;
            pos_a_end[hd] = link_length;
        }
        else{ 
            /*Move the motor to the next link on the filament
             *At the projected new position along that filament*/
            l_index[hd] = l_index[hd] - 1;
            pos_a_end[hd] = pos - link_length;
        }
    }
    else if (pos < 0) { //this shouldn't be possible if vm > 0
        if (l_index[hd] == (actin_network->get_filament(f_index[hd])->get_nlinks() - 1)){ // the pointed end of the filament
            pos_a_end[hd]=0; //move head to pointed end
        }
        else{ 
            /*Move the motor to the previous link on the filament
             *At the projected new position along that filament*/
            l_index[hd] = l_index[hd] + 1;
            pos_a_end[hd] = pos + actin_network->get_llength(f_index[hd],l_index[hd]);    
        }
    }   
    else {
        pos_a_end[hd] = pos;
    }
       
}


void motor::update_position_attached(int hd){

    double posx = actin_network->get_end(f_index[hd],l_index[hd])[0]-pos_a_end[hd]*actin_network->get_direction(f_index[hd],l_index[hd])[0];
    double posy = actin_network->get_end(f_index[hd],l_index[hd])[1]-pos_a_end[hd]*actin_network->get_direction(f_index[hd],l_index[hd])[1];

    array<double, 2> newpos = boundary_check(hd, posx, posy);
    
    hx[hd] = newpos[0];
    hy[hd] = newpos[1];

}

// Using the lever rule to propagate force as outlined in Nedelec F 2002

void motor::actin_update_hd(int hd, array<double, 2> f)
{
    double pos_ratio = pos_a_end[hd]/actin_network->get_llength(f_index[hd], l_index[hd]);
    actin_network->update_forces(f_index[hd], l_index[hd],   f[0] *    pos_ratio , f[1] *    pos_ratio );
    actin_network->update_forces(f_index[hd], l_index[hd]+1, f[0] * (1-pos_ratio), f[1] * (1-pos_ratio));
}


void motor::actin_update()
{
    if (state[0]==1 && state[1]==0) this->actin_update_hd(0, force);
    if (state[0]==0 && state[1]==1) this->actin_update_hd(1, {-force[0], -force[1]});
    if (state[0]==1 && state[1] == 1)
    {
        //cout <<"angular update"<<endl;
        this->actin_update_hd(0, force);
        this->actin_update_hd(1, {-force[0], -force[1]});
        this->actin_update_angular();
    }
}

void motor::actin_update_angular()
{
      array<double, 2> delr1, delr2;
      double f1[2], f3[2], f2[2], f4[2];
      double rsq1,rsq2,r1,r2,c,s,a,a11,a12,a22;
      delr1 = actin_network->get_filament(f_index[0])->get_link(l_index[0])->get_neg_disp();
      delr2 = actin_network->get_filament(f_index[1])->get_link(l_index[1])->get_neg_disp();
//      double fa[2],fb[2],fc[2],fd[2];
      // 1st bond
      //delr1 = links[0]->get_neg_disp();
      rsq1  = delr1[0]*delr1[0] + delr1[1]*delr1[1];
      r1    = sqrt(rsq1);
      
      // 2nd bond
      //delr2 = links[n+1]->get_disp();
      rsq2  = delr2[0]*delr2[0] + delr2[1]*delr2[1];
      r2    = sqrt(rsq2);
        
//      double actin_0_x = actin_network->get_filament(f_index[0])->get_actin(l_index[0])->get_xcm();
//      double actin_0_y = actin_network->get_filament(f_index[0])->get_actin(l_index[0])->get_ycm();
//      double actin_1_x = actin_network->get_filament(f_index[1])->get_actin(l_index[1])->get_xcm();
//      double actin_1_y = actin_network->get_filament(f_index[1])->get_actin(l_index[1])->get_ycm();
//    double b_head_x_0 = actin_0_x - hx[0];
//    double b_head_y_0 = actin_0_y - hy[0];
//    double b_head_x_1 = actin_1_x - hx[1];
//    double b_head_y_1 = actin_1_y - hy[1];
    
//      double pos_ratio_0 = pos_a_end[0]/actin_network->get_llength(f_index[0], l_index[0]);
//      double pos_ratio_1 = pos_a_end[1]/actin_network->get_llength(f_index[1], l_index[1]);
      
    //cout<<"\npos_ratio_0: "<<pos_ratio_0<<" - pos_ratio_1: "<<pos_ratio_1<<endl;
      // angle (cos and sin)
      c = (delr1[0]*delr2[0] + delr1[1]*delr2[1]) / (r1*r2);
      
      if (c > 1.0) c = 1.0;
      if (c < -1.0) c = -1.0;

      s = sqrt(1.0 - c*c);
      if (s < maxSmallAngle) s = maxSmallAngle;

      // force
      //cout<<"\nmax_small_angle: "<<a;
      a   = -0.011 * (acos(c) - 0) / s; //Note, acos(c) in radians
      a11 = a*c / rsq1;
      a12 = a / (r1*r2);
      a22 = a*c / rsq2;
    
    //cout<<"\na, a11, a22: "<<a<<" "<<a11<<" "<<a22;
    
      f1[0] = a11*delr1[0] - a12*delr2[0];
      f1[1] = a11*delr1[1] - a12*delr2[1];
      f2[0] = -a11*delr1[0] + a12*delr2[0];
      f2[1] = -a11*delr1[1] + a12*delr2[1];
      f3[0] = a22*delr2[0] - a12*delr1[0];
      f3[1] = a22*delr2[1] - a12*delr1[1];
      f4[0] = -a22*delr2[0] + a12*delr1[0];
      f4[1] = -a22*delr2[1] + a12*delr1[1];
    
//    double force_a = sqrt(f1[0]*f1[0] + f1[1]*f1[1]);
//    double force_b = sqrt(f2[0]*f2[0] + f2[1]*f2[1]);
//    double force_c = sqrt(f3[0]*f3[0] + f3[1]*f3[1]);
//    double force_d = sqrt(f4[0]*f4[0] + f4[1]*f4[1]);
//    
//    fa[0] = (force_a + force_b) * (pos_ratio_0)*f1[0]/force_a;
//    fa[1] = (force_a + force_b) * (pos_ratio_0)*f1[1]/force_a;
//    fb[0] = (force_a + force_b)*f2[0] * (1-pos_ratio_0)/force_b;
//    fb[1] = (force_a + force_b)*f2[1] * (1-pos_ratio_0)/force_b;
//    fc[0] = (force_c + force_d) * (pos_ratio_1)*f3[0]/force_c;
//    fc[1] = (force_c + force_d) * (pos_ratio_1)*f3[1]/force_c;
//    fd[0] = (force_c + force_d)*f4[0] * (1-pos_ratio_1)/force_d;
//    fd[1] = (force_c + force_d)*f4[1] * (1-pos_ratio_1)/force_d;
    
    //cout<<"\nforce_a: "<<force_a<<endl;
    //cout<<"\nDEBUG: f1x, f1y, f2x, f2y, f3x, f3y, f4x, f4y = "<<f1[0]<<" , "<<f1[1]<<" , "<<f2[0]<<" , "<<f2[1]<<" , "<<f3[0]<<" , "<<f3[1]<<" , "<<f4[0]<<" , "<<f4[1];
    //cout<<"\nDEBUG: f1x, f1y, f2x, f2y, f3x, f3y, f4x, f4y = "<<fa[0]<<" , "<<fa[1]<<" , "<<fb[0]<<" , "<<fb[1]<<" , "<<fc[0]<<" , "<<fc[1]<<" , "<<fd[0]<<" , "<<fd[1];
    

    
    actin_network->update_forces(f_index[0], l_index[0], f1[0], f1[1]);
    actin_network->update_forces(f_index[0], l_index[0]+1, f2[0], f2[1]);
    actin_network->update_forces(f_index[1], l_index[1], f3[0], f3[1]);
    actin_network->update_forces(f_index[1], l_index[1]+1, f4[0], f4[1]);
    
 //   actin_network->update_forces(f_index[0], l_index[0], fa[0], fa[1]);
 //   actin_network->update_forces(f_index[0], l_index[0]+1, fb[0], fb[1]);
 //   actin_network->update_forces(f_index[1], l_index[1], fc[0], fc[1]);
 //   actin_network->update_forces(f_index[1], l_index[1]+1, fd[0], fd[1]);
    //cout<<"\nDEBUG: f_index0, l_index0, f_index1, l_index1 = "<<f_index[0]<<" , "<<l_index[0]<<" , "<<f_index[1]<<" , "<<l_index[1];

}

void motor::detach_head(int hd, array<double, 2> newpos)
{
   
    state[hd]=0;
    f_index[hd]=-1;
    l_index[hd]=-1;
    pos_a_end[hd]=0;
    
    hx[hd] = newpos[0];
    hy[hd] = newpos[1];

}

void motor::detach_head_without_moving(int hd)
{
   
    state[hd]=0;
    f_index[hd]=-1;
    l_index[hd]=-1;
    pos_a_end[hd]=0;
    
}


array<int, 2> motor::get_f_index(){
    return f_index;
}


array<int, 2> motor::get_l_index(){
    return l_index;
}


array<double, 2> motor::get_pos_a_end(){
    return pos_a_end;
}


array<double, 2> motor::get_force(){
    return force;
}


double motor::get_stretching_energy(){
    return (force[0]*force[0]+force[1]*force[1])/(2*mk);
}


double motor::get_stretching_energy_fene()
{
    double ext = abs(mld - hypot(disp[0], disp[1]));
    
    if (max_ext - ext > eps_ext )
        return -0.5*mk*max_ext*max_ext*log(1-(ext/max_ext)*(ext/max_ext));
    else
        return 0.25*mk*ext*ext*(max_ext/eps_ext);
    
}


double motor::get_kinetic_energy(){
    return kinetic_energy;
}


string motor::to_string()
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
}


string motor::write()
{
    return std::to_string(hx[0]) + "\t" + std::to_string(hy[0])
        +  "\t" + std::to_string(disp[0]) + "\t" + std::to_string(disp[1]) 
        +  "\t" + std::to_string(f_index[0]) + "\t" + std::to_string(f_index[1]) 
        +  "\t" + std::to_string(l_index[0]) + "\t" + std::to_string(l_index[1]) + "\n";
}
