/*
    Copyright (C) <2014-2020>  <PingzhouMing>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */
//
// Class  : NeutronTransport v1.0 (~750 lines of c++ code including comments)
// Date   : June 2009
//
// Simulates a Neutron Transport System for a spherical bulk (medium) containing U235 and U238 nuclei.
//
// Thermal neutrons are relased from the center of a shperical bulk of radius m_radius and mass m_mass
// containing only two nuclei U235 and U238. Using Monte Carlo methods, program simulates
// the trajectories (i.e. the positions) of the secondary neutrons propagating in the bulk.
// The secondaries can be prompt neutrons (from fission reaction) or scattered neutrons.
//
// At the end, end() method outputs the average neutron multiplication factor (<keff>) of the
// system defined as the ratio:
//
//   keff = (# neutrons generated in one generation) / (# neutrons generated in preceding generation)
//
//   keff < 1 ==> results in sub-critial mass   (means no chain reaction occurs)
//   keff = 1 ==> results in critial mass       (means a chain reaction can start)
//   keff > 1 ==> results in super-critial mass (means a chain reaction definitely starts)
//
// You can find more information at: http://www1.gantep.edu.tr/~bingul/simulation/fission/
// and a main (driver) program at:   http://www1.gantep.edu.tr/~bingul/simulation/fission/fission.cc
//
// Author : Dr. Ahmet Bing <bingul(at)gantep.edu.tr>
// Modified by PingzhouMing@2014 <mingpz@qq.com>
//          2014.11.23,   add memset() to ensure init conditions!
//          2014.12.15,   enable 3D neutron simulation
//          2015.05.02,   simulate nine cube assemblies to simplified model
//          2015.06.20,   use DMAT library to optimize simulation
//
#ifndef FISSION_H_
#define FISSION_H_

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <string>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <vector>
#include <fstream>
#include "unistd.h"
#include "../src/Matrix_hpc.h"

using namespace std;

const int ZONE = 9;
typedef struct Neutron {
    double x;
    double y;
    double z;
    double ux;
    double uy;
    double uz;
    double e;
}Neutron_t;

class NeutronTransport3D
{
private:
    // *** private variables and methods starts with m_  ***
    double m_PI, m_A0, m_barn, m_A235, m_A238, m_N235[ZONE], m_N238[ZONE],
       m_density[ZONE], m_purity[ZONE], m_eAvr, m_eCharge,
       m_micro235[ZONE][5], m_macro235[ZONE][5], m_micro238[ZONE][5], m_macro238[ZONE][5],
       m_macroAll[ZONE][5], m_P[5], m_lambda[ZONE], m_heat, m_En, m_radius, m_mass[ZONE],
       m_s235[5][20], m_s238[5][20], m_Energy[20], m_eThermal,
       m_pAbs[ZONE], m_pCap[ZONE], m_pFis[ZONE], m_p235[ZONE], m_pSct[ZONE], m_pEla[ZONE], m_pInE[ZONE];
    bool   m_showDetails;
    int    m_evtMax, m_evtNum, m_genNum, m_genMax, m_nMax, **m_nNeutron[ZONE], m_hit, *m_Count[ZONE];

    int number_of_enter[ZONE];
    // To save position of neutrons
    int p_start[ZONE], p_enter[ZONE];
    double *px[ZONE], *py[ZONE], *pz[ZONE], *pe[ZONE], *pux[ZONE], *puy[ZONE], *puz[ZONE];
    Matrix_hpc<double> px_opt;
    Matrix_hpc<double> py_opt;
    Matrix_hpc<double> pz_opt;
    Matrix_hpc<double> pe_opt;
    Matrix_hpc<double> pux_opt;
    Matrix_hpc<double> puy_opt;
    Matrix_hpc<double> puz_opt;

    int    m_randomSeed;
    double m_randomFlat(void);
    int    m_randomNubar(double, bool);
    double m_randomPathLength(int);
    double m_randomPromptEnergy(void);
    double m_randomScatteredEnergy(double, double, bool);

    int    m_chain_new(int section, double x, double y, double z, double e, double ux, double uy, double uz);
    int    m_chain_new_opt(int section, double x, double y, double z, double e, double ux, double uy, double uz);
    int    m_nint(double);
    void   m_getData(void);
    void   m_getProbabilities(double);
    bool   m_isOutside(double, double, double);
    // Velocity keeps the same coordinates
    int    m_enter_other_section(int source, double&, double&, double&);

    time_t m_timeStart, m_timeEnd;

public:
    NeutronTransport3D(int, int, double pure[ZONE], bool, int);

    void   start();
    void   run(int type);
    void   execute_opt(int k);
    void   execute(int k);
    void   end();
};

NeutronTransport3D::NeutronTransport3D(int max_evt, int max_gen, double pur[ZONE], bool detail, int seed)
{
  //*** Parameters used in the program ***              // Description                                Unit
                                                        // ----------------------------------------   ----------
  m_PI      = 3.14159265358979323846;                   // Number pi                                  -
  m_A0      = 6.02214179e+23;                           // Avagadros number                           atoms/mole
  m_barn    = 1.0e-24;                                  // A conversion factor between barn and cm2   cm2
  m_A235    = 235.0439299;                              // Atomic weigth of U235                      u
  m_A238    = 238.0507826;                              // Atomic weigth of U235                      u
  for (int i=0; i<ZONE; i++) {
      m_density[i] = 18.75;                                    // mass density of Uranium                    g/cm3
      if(pur[i] < 0.9) m_density[i] = 19.05;                   // mass density for purity < 90%              g/cm3
      m_purity[i]  = pur[i];                                      // purity of the U235                         -
      m_N235[i]    = m_purity[i]*m_density[i]*m_A0/m_A235;           // Atom density of U235                       1/cm3
      m_N238[i]    = (1.0-m_purity[i])*m_density[i]*m_A0/m_A238;     // Atom density of U238
      m_mass[i]    = m_density[i] * 20*20*100;                    // mass of the cube
  }
  m_eCharge = 1.602e-19;                                // electronic charge                          C
  m_eAvr    = 2.0e+8 * m_eCharge;                       // Average energy released per fission        J1/cm3
  m_heat    = 0.0;                                      // Average kinetic energy per fission         J                    g
  m_eThermal= 0.0253;                                   // Kinetic energy of thermal neutrons         eV

  //*** Neutron counting stuff ***
  m_hit    = 0;           // see m_printDetails() routine
  m_nMax   = 200000;       // Maximum number of neutrons in one generation
  m_evtNum = 0;           // Event Number counter (upto m_evtMax)
  m_genNum = 0;           // Generation number counter (upto m_genMax)
  m_evtMax = max_evt;     // Maximum number of events (default is 1)
  m_genMax = max_gen;     // Maximum number of generations allowed (default is 20)
  m_showDetails = detail; // To show details of the each event if detail = true (default is false)

  //*** total number of neutrons generated ***
#if 0
  for (int i=0; i<ZONE; i++) {
      m_nNeutron[i] = new int* [m_evtMax];
  }
  for(int i=0; i<m_evtMax; i++) {
       for(int j=0; j<ZONE; j++) {
           *(m_nNeutron[j]+i) = new int [m_genMax];
           memset(*(m_nNeutron[j]+i), 0, sizeof(int)*m_genMax);
       }
  } //
#endif
  for (int i=0; i<ZONE; i++) {
      m_Count[i] = new int [m_genMax];
      m_Count[i][0] = max_evt;
      for (int j=1; j<m_genMax; j++) {
          m_Count[i][j] = 0;
      }
    }

  //*** initilize the random number generator ***
  if( !seed ) m_randomSeed = (int) time(NULL);
  else        m_randomSeed = seed;

  //*** to calculate program execution time count number of seconds ***
  time(&m_timeStart);
}

void NeutronTransport3D::run(int type)
{
    if (type == 0) {
        for (int k=0; k<m_genMax-1; k++) {
            execute(k);
        }
    } else {
        for (int k=0; k<m_genMax-1; k++) {
            execute_opt(k);
        }
    }

    if (type == 0) {
        for (int i=0; i<ZONE; i++) {
            if (px[i] != NULL)      delete [] px[i];
            if (py[i] != NULL)      delete [] py[i];
            if (pz[i] != NULL)      delete [] pz[i];
            if (pe[i] != NULL)      delete [] pe[i];
            if (pux[i] != NULL)     delete [] pux[i];
            if (puy[i] != NULL)     delete [] puy[i];
            if (puz[i] != NULL)     delete [] puz[i];
        }
    }
}

void NeutronTransport3D::end()
{
    double  keff[m_genMax], kerr[m_genMax];

    cout << endl << "--- Neutron Transport stoped. --------------------------------------------" << endl;
    cout << "Total Number of Events: " << m_evtMax << endl;
    cout << "Total Generation      : " << m_genMax << endl << endl;
    for (int i=0; i<ZONE; i++) {
        cout<<"Assemble "<< i+1 <<":"<<endl;
#if 0
        cout << "Generation #    nold          nnew              keff          kerr   " << endl;
        cout << "============    ==========    ==========        =======       =======" << endl;
#endif
        for(int j=1; j<m_genMax; j++)
        {
          double ratio = double(m_Count[i][j]) / m_Count[i][j-1];

          if(ratio<100.)
          {
               keff[j] = ratio;
               kerr[j] = sqrt( 1.0/m_Count[i][j] + 1.0/m_Count[i][j-1]) * ratio;
          }
#if 0
          cout << setw(12) << j << setw(14) << m_Count[i][j-1] << setw(14) << m_Count[i][j]
               << fixed << setprecision(5) << setw(15) << keff[j]
               << fixed << setprecision(5) << setw(14) << kerr[j]
               << endl;
#endif
        }

        double sum = 0, ne = 0;
        for(int i=0; i<m_genMax; i++){
           if(keff[i]) {sum += keff[i]; ne++;}
        }
        double keffavr = sum/ne;

        cout << "<keff> = " << fixed << setprecision(5)<<keffavr<<endl<<endl;
    }
}

int NeutronTransport3D::m_enter_other_section(int source,
        double& x, double& y, double& z)
{
    if (z >= 50 || z <= -50) {
        return -1;
    }
    if (source == 0) {
        if (x>=10 && (y<10 && y>-10)) {
            x = x-20;
            return 1;
        }
        if (x>=10 && y>=10) {
            x = x-20;
            y = y-20;
            return 2;
        }
        if (x>=10 && y<=-10) {
            x = x-20;
            y = y+20;
            return 8;
        }
        if (x<=-10 && (y<10 && y>-10)) {
            x = x+20;
            return 5;
        }
        if (x<=-10 && y>=10) {
            x = x+20;
            y = y-20;
            return 4;
        }
        if (x<=-10 && y<=-10) {
            x = x+20;
            y = y+20;
            return 6;
        }
        if (y>=10 && (x<10 && x>-10) ) {
            y = y-20;
            return 3;
        }
        if (y>=10 && x>=10) {
            x = x-20;
            y = y-20;
            return 2;
        }
        if (y>=10 && x<=-10) {
            x = x+20;
            y = y-20;
            return 4;
        }
        if (y<-10 && (x<10 && x>-10)) {
            y = y+20;
            return 7;
        }
        if (y<-10  && x>=10 ) {
            x = x-20;
            y = y+20;
            return 8;
        }
        if (y<-10 && x<=-10) {
            x = x+20;
            y = y+20;
            return 6;
        }
    }
    if (source == 1) {
        if (x>=10) {
            return -1;
        }
        if (x<=-10 && (y<10 && y>-10)) {
            x = x+20;
            return 0;
        }
        if (x<=-10 && y>=10) {
            x = x+20;
            y = y-20;
            return 3;
        }
        if (x<=-10 && y<=-10) {
            x = x+20;
            y = y+20;
            return 7;
        }
        if (y>=10 && (x<10 && x>-10) ) {
            y = y-20;
            return 2;
        }
        if (y>=10 && x>=10) {
            return -1;
        }
        if (y>=10 && x<=-10) {
            x = x+20;
            y = y-20;
            return 3;
        }
        if (y<-10 && (x<10 && x>-10)) {
            y = y+20;
            return 8;
        }
        if (y<-10  && x>=10 ) {
            return -1;
        }
        if (y<-10 && x<=-10) {
            x = x+20;
            y = y+20;
            return 7;
        }
    }
    if (source == 2) {
        if (x<=-10 && (y<10 && y>-10)) {
            x = x+20;
            return 3;
        }
        if (x<=-10 && y>=10) {
            return -1;
        }
        if (x<=-10 && y<=-10) {
            x = x+20;
            y = y+20;
            return 0;
        }
        if (y<-10 && (x<10 && x>-10)) {
            y = y+20;
            return 1;
        }
        if (y<-10  && x>=10 ) {
            return -1;
        }
        if (y<-10 && x<=-10) {
            x = x+20;
            y = y+20;
            return 0;
        }
    }
    if (source == 3) {
        if (x>=10 && (y<10 && y>-10)) {
            x = x-20;
            return 2;
        }
        if (x>=10 && y>=10) {
            return -1;
        }
        if (x>=10 && y<=-10) {
            x = x-20;
            y = y+20;
            return 1;
        }
        if (x<=-10 && (y<10 && y>-10)) {
            x = x+20;
            return 4;
        }
        if (x<=-10 && y>=10) {
            return -1;
        }
        if (x<=-10 && y<=-10) {
            x = x+20;
            y = y+20;
            return 5;
        }
        if (y<-10 && (x<10 && x>-10)) {
            y = y+20;
            return 0;
        }
        if (y<-10  && x>=10 ) {
            x = x-20;
            y = y+20;
            return 1;
        }
        if (y<-10 && x<=-10) {
            x = x+20;
            y = y+20;
            return 5;
        }
    }
    if (source == 4) {
        if (x>=10 && (y<10 && y>-10)) {
            x = x-20;
            return 3;
        }
        if (x>=10 && y>=10) {
            return -1;
        }
        if (x>=10 && y<=-10) {
            x = x-20;
            y = y+20;
            return 0;
        }
        if (y<-10 && (x<10 && x>-10)) {
            y = y+20;
            return 5;
        }
        if (y<-10  && x>=10 ) {
            x = x-20;
            y = y+20;
            return 0;
        }
        if (y<-10 && x<=-10) {
            return -1;
        }
    }
    if (source == 5) {
        if (x>=10 && (y<10 && y>-10)) {
            x = x-20;
            return 0;
        }
        if (x>=10 && y>=10) {
            x = x-20;
            y = y-20;
            return 3;
        }
        if (x>=10 && y<=-10) {
            x = x-20;
            y = y+20;
            return 7;
        }
        if (y>=10 && (x<10 && x>-10) ) {
            y = y-20;
            return 4;
        }
        if (y>=10 && x>=10) {
            x = x-20;
            y = y-20;
            return 3;
        }
        if (y>=10 && x<=-10) {
            return -1;
        }
        if (y<-10 && (x<10 && x>-10)) {
            y = y+20;
            return 6;
        }
        if (y<-10  && x>=10 ) {
            x = x-20;
            y = y+20;
            return 7;
        }
        if (y<-10 && x<=-10) {
            return -1;
        }
    }
    if (source == 6) {
        if (x>=10 && (y<10 && y>-10)) {
            x = x-20;
            return 7;
        }
        if (x>=10 && y>=10) {
            x = x-20;
            y = y-20;
            return 0;
        }
        if (x>=10 && y<=-10) {
            return -1;
        }
        if (y>=10 && (x<10 && x>-10) ) {
            y = y-20;
            return 5;
        }
        if (y>=10 && x>=10) {
            x = x-20;
            y = y-20;
            return 0;
        }
        if (y>=10 && x<=-10) {
            return -1;
        }
    }
    if (source == 7) {
        if (x>=10 && (y<10 && y>-10)) {
            x = x-20;
            return 8;
        }
        if (x>=10 && y>=10) {
            x = x-20;
            y = y-20;
            return 1;
        }
        if (x>=10 && y<=-10) {
            return -1;
        }
        if (x<=-10 && (y<10 && y>-10)) {
            x = x+20;
            return 6;
        }
        if (x<=-10 && y>=10) {
            x = x+20;
            y = y-20;
            return 5;
        }
        if (x<=-10 && y<=-10) {
            return -1;
        }
        if (y>=10 && (x<10 && x>-10) ) {
            y = y-20;
            return 0;
        }
        if (y>=10 && x>=10) {
            x = x-20;
            y = y-20;
            return 1;
        }
        if (y>=10 && x<=-10) {
            x = x+20;
            y = y-20;
            return 5;
        }
    }
    if (source == 8) {
        if (x<=-10 && (y<10 && y>-10)) {
            x = x+20;
            return 7;
        }
        if (x<=-10 && y>=10) {
            x = x+20;
            y = y-20;
            return 0;
        }
        if (x<=-10 && y<=-10) {
            return -1;
        }
        if (y>=10 && (x<10 && x>-10) ) {
            y = y-20;
            return 1;
        }
        if (y>=10 && x>=10) {
            return -1;
        }
        if (y>=10 && x<=-10) {
            x = x+20;
            y = y-20;
            return 0;
        }
    }

    return -1;
}

void NeutronTransport3D::execute_opt(int k)
{
    memset(number_of_enter, 0, sizeof(int)*ZONE);
    if (k == 0) {
        for (int i=0; i<ZONE; i++) {
              p_start[i] = 0;  p_enter[i] = m_evtMax;
              px_opt = px_opt.newsize(ZONE, m_nMax);
              py_opt = py_opt.newsize(ZONE, m_nMax);
              pz_opt = pz_opt.newsize(ZONE, m_nMax);
              pe_opt = pe_opt.newsize(ZONE, m_nMax);
              pux_opt = pux_opt.newsize(ZONE, m_nMax);
              puy_opt = puy_opt.newsize(ZONE, m_nMax);
              puz_opt = puz_opt.newsize(ZONE, m_nMax);
          // initial positions of the neutrons for the first generation
          for(int j=0; j<m_evtMax; )
          {
             double phi = 2.0*m_PI*m_randomFlat();
             double costh = 2.0*m_randomFlat()-1.0;
             double sinth = sqrt(1.0-costh*costh);
             double rad = 10 * sqrt(27) * m_randomFlat();  // cube组件采用球坐标系
             px_opt(i,j) = rad*sinth*cos(phi);
             py_opt(i,j) = rad*sinth*sin(phi);
             pz_opt(i,j) = rad*costh;
             if ( m_isOutside(px_opt(i,j), py_opt(i,j), pz_opt(i,j)) ) {
                 continue;
             }
             pux_opt(i,j) = -px_opt(i,j);
             puy_opt(i,j) = -py_opt(i,j);
             puz_opt(i,j) = -pz_opt(i,j);

             //x[j] = y[j] = z[j] = ux[j] = uy[j] = uz[j] = 0.0;
             pe_opt(i,j) = m_eThermal;
             j++;
              }
          }
      }

      // event loop **********************************************************
      for (int i=0; i<ZONE; i++) {
          int n = m_Count[i][k];
          for(int m = 0; m<n; m++) {
            // Each event starts with n = 1 neutron

            // Total # of neutrons generated for this event
            // generation loop
            m_Count[i][k+1] += m_chain_new_opt(i, px_opt(i, p_start[i]), py_opt(i, p_start[i]),
                    pz_opt(i, p_start[i]), pe_opt(i, p_start[i]), pux_opt(i, p_start[i]),
                    puy_opt(i, p_start[i]), puz_opt(i, p_start[i]));;  // k+1 generation
            p_start[i]++;
            if (p_start[i] >= m_nMax) {
                p_start[i] %= m_nMax;
            }
          } // end of event loop *************************************************
          sleep(1);
          printf("%g\n", m_randomFlat());
      }

      for (int i=0; i<ZONE; i++) {
          m_Count[i][k+1] += number_of_enter[i];
          // cout<<"n = "<<n<<" p_start = "<<p_start[i]<<" p_enter = "<<p_enter[i]<<endl;
#if 1
          cout << "Assemble: " << setw(7) << i+1
               << " | Generation: " << setw(5) << k+1
               << " | Total # of genenerated neutrons: " << m_Count[i][k+1]
               << endl;
#endif
      }
}

void NeutronTransport3D::execute(int k)
{
    memset(number_of_enter, 0, sizeof(int)*ZONE);
    if (k == 0) {
        for (int i=0; i<ZONE; i++) {
              p_start[i] = 0;  p_enter[i] = m_evtMax;
              px[i] = new double  [m_nMax];
              py[i] = new double  [m_nMax];
              pz[i] = new double  [m_nMax];
              pe[i] = new double  [m_nMax];
              pux[i] = new double [m_nMax];
              puy[i] = new double [m_nMax];
              puz[i] = new double [m_nMax];
          // initial positions of the neutrons for the first generation
          for(int j=0; j<m_evtMax; )
          {
             double phi = 2.0*m_PI*m_randomFlat();
             double costh = 2.0*m_randomFlat()-1.0;
             double sinth = sqrt(1.0-costh*costh);
             double rad = 10 * sqrt(27) * m_randomFlat();  // cube组件采用球坐标系
             px[i][j] = rad*sinth*cos(phi);
             py[i][j] = rad*sinth*sin(phi);
             pz[i][j] = rad*costh;
             if ( m_isOutside(px[i][j], py[i][j], pz[i][j]) ) {
                 continue;
             }
             pux[i][j] = -px[i][j];
             puy[i][j] = -py[i][j];
             puz[i][j] = -pz[i][j];

             //x[j] = y[j] = z[j] = ux[j] = uy[j] = uz[j] = 0.0;
             pe[i][j] = m_eThermal;
             j++;
              }
          }
      }

      // event loop **********************************************************
      for (int i=0; i<ZONE; i++) {
          int n = m_Count[i][k];
          for(int m = 0; m<n; m++) {
            // Each event starts with n = 1 neutron

            // Total # of neutrons generated for this event
            // generation loop
            m_Count[i][k+1] += m_chain_new(i, px[i][p_start[i]], py[i][p_start[i]], pz[i][p_start[i]], pe[i][p_start[i]], pux[i][p_start[i]], puy[i][p_start[i]], puz[i][p_start[i]]);;  // k+1 generation
            p_start[i]++;
            if (p_start[i] >= m_nMax) {
                p_start[i] %= m_nMax;
            }
          } // end of event loop *************************************************
          sleep(1);
          printf("%g\n", m_randomFlat());
      }

      for (int i=0; i<ZONE; i++) {
          m_Count[i][k+1] += number_of_enter[i];
          // cout<<"n = "<<n<<" p_start = "<<p_start[i]<<" p_enter = "<<p_enter[i]<<endl;
#if 1
          cout << "Assemble: " << setw(7) << i+1
               << " | Generation: " << setw(5) << k+1
               << " | Total # of genenerated neutrons: " << m_Count[i][k+1]
               << endl;
#endif
      }
}

int NeutronTransport3D::m_chain_new_opt
(int section, double x, double y, double z, double e, double ux, double uy, double uz)
{
    bool isAbsorbed  = false;
    bool isScattered = false;
    bool isFission   = false;
    bool isCaptured  = false;
    bool isElastic   = false;
    bool isInElastic = false;
    bool isPossible  = true;

    m_getProbabilities( e );
    int ng = 0;
    if( m_randomFlat() < m_pAbs[section] ){
        isAbsorbed = true;
        if( m_randomFlat() < m_pFis[section] ) isFission  = true; // fission
        else                          isCaptured = true; // capture
    } // if scattering occur
    else{
        isScattered = true;
        if( m_randomFlat() < m_pEla[section] ) isElastic   = true; // elastic scatt.
        else                          isInElastic = true; // in-elastic scatt.
    } //************************************************************************

    if(isFission) {
        // is the nucleus U235 ?
        bool isU235 = false;
        if(m_randomFlat() < m_p235[section]) isU235 = true;

        // Get number of prompt neutrons produced for the fission
        int np = m_randomNubar(e, isU235);

        for(int l=1; l<=np; l++) {
         // Get a random free path length
         double d = m_randomPathLength(section);

         // Send the prompt neutron to a isotropiaclly random direction in space
         double phi  = 2.0*m_PI*m_randomFlat();
         double theta= acos(2.0*m_randomFlat()-1.0);
         double uux = sin(theta)*cos(phi);
         double uuy = sin(theta)*sin(phi);
         double uuz = cos(theta);
         double xx  = x + d*uux;
         double yy  = y + d*uuy;
         double zz  = z + d*uuz;

         // Assign a random kinetic energy to the prompt neutron
         double ee = m_randomPromptEnergy();

         // Check if the neutron is outside of the bulk
         if( m_isOutside(xx, yy, zz) ) {
             int index = m_enter_other_section(section, xx, yy, zz);
             if (index >= 0) {  // >=0，因为组件数目在计算机中是从0开始计数
                 px_opt(index, p_enter[index]) = xx;
                 py_opt(index, p_enter[index]) = yy;
                 pz_opt(index, p_enter[index]) = zz;
                 pe_opt(index, p_enter[index]) = ee;
                 pux_opt(index, p_enter[index]) = uux;
                 puy_opt(index, p_enter[index]) = uuy;
                 puz_opt(index, p_enter[index]) = uuz;
                 p_enter[index]++;
                 if (p_enter[index] >= m_nMax) {
                     p_enter[index] %= m_nMax;
                 }
                 number_of_enter[index]++;
             }
         } else{
                 px_opt(section, p_enter[section]) = xx;
                 py_opt(section, p_enter[section]) = yy;
                 pz_opt(section, p_enter[section]) = zz;
                 pe_opt(section, p_enter[section]) = ee;
                 pux_opt(section, p_enter[section]) = uux;
                 puy_opt(section, p_enter[section]) = uuy;
                 puz_opt(section, p_enter[section]) = uuz;
                 p_enter[section]++;
                 if (p_enter[section] >= m_nMax) {
                     p_enter[section] %= m_nMax;
                 }
                 ng++;
             }
         } // end of prompt neutron loop
      } // end of isFission
    else if(isScattered) {
        // Get a random free path length for this energy
          double d = m_randomPathLength(section);

          // Send the scattered neutron to a isotropiaclly random direction in space
          double phi  = 2.0*m_PI*m_randomFlat();
          double theta= acos(2.0*m_randomFlat()-1.0);
          double uux = sin(theta)*cos(phi);
          double uuy = sin(theta)*sin(phi);
          double uuz = cos(theta);
          double xx  = x + d*uux;
          double yy  = y + d*uuy;
          double zz  = z + d*uuz;

          // Scattering angle
          // i.e. opening angle between new and old direction (it is found from dot product)
          double cosThetaS = ux*uux + uy*uuy + uz*uuz;

          // Determine the random kinetic energy to this scattered neutron
          double ee = m_randomScatteredEnergy(e, cosThetaS, isElastic);
          if(ee<0.0) isPossible = false;

          // Check if the neutron is outside of the bulk
          if( m_isOutside(xx, yy, zz) ) {
              int index = m_enter_other_section(section, xx, yy, zz);
              if (index >= 0) {
                  px_opt(section, p_enter[section]) = xx;
                   py_opt(section, p_enter[section]) = yy;
                   pz_opt(section, p_enter[section]) = zz;
                   pe_opt(section, p_enter[section]) = ee;
                   pux_opt(section, p_enter[section]) = uux;
                   puy_opt(section, p_enter[section]) = uuy;
                   puz_opt(section, p_enter[section]) = uuz;
                  p_enter[index]++;
                  if (p_enter[index] >= m_nMax) {
                      p_enter[index] %= m_nMax;
                  }
                  number_of_enter[index]++;
              }
          } else if(isPossible){
              px_opt(section, p_enter[section]) = xx;
               py_opt(section, p_enter[section]) = yy;
               pz_opt(section, p_enter[section]) = zz;
               pe_opt(section, p_enter[section]) = ee;
               pux_opt(section, p_enter[section]) = uux;
               puy_opt(section, p_enter[section]) = uuy;
               puz_opt(section, p_enter[section]) = uuz;
              p_enter[section]++;
              if (p_enter[section] >= m_nMax) {
                  p_enter[section] %= m_nMax;
              }
               ng++;
          }
          else {
              (void)0;
          }
       } else if(isCaptured)
       {
             (void)0;
       }
    return ng;
}

int NeutronTransport3D::m_chain_new
(int section, double x, double y, double z, double e, double ux, double uy, double uz)
{
    bool isAbsorbed  = false;
    bool isScattered = false;
    bool isFission   = false;
    bool isCaptured  = false;
    bool isElastic   = false;
    bool isInElastic = false;
    bool isPossible  = true;

    m_getProbabilities( e );
    int ng = 0;
    if( m_randomFlat() < m_pAbs[section] ){
        isAbsorbed = true;
        if( m_randomFlat() < m_pFis[section] ) isFission  = true; // fission
        else                          isCaptured = true; // capture
    } // if scattering occur
    else{
        isScattered = true;
        if( m_randomFlat() < m_pEla[section] ) isElastic   = true; // elastic scatt.
        else                          isInElastic = true; // in-elastic scatt.
    } //************************************************************************

    if(isFission) {
        // is the nucleus U235 ?
        bool isU235 = false;
        if(m_randomFlat() < m_p235[section]) isU235 = true;

        // Get number of prompt neutrons produced for the fission
        int np = m_randomNubar(e, isU235);

        for(int l=1; l<=np; l++) {
         // Get a random free path length
         double d = m_randomPathLength(section);

         // Send the prompt neutron to a isotropiaclly random direction in space
         double phi  = 2.0*m_PI*m_randomFlat();
         double theta= acos(2.0*m_randomFlat()-1.0);
         double uux = sin(theta)*cos(phi);
         double uuy = sin(theta)*sin(phi);
         double uuz = cos(theta);
         double xx  = x + d*uux;
         double yy  = y + d*uuy;
         double zz  = z + d*uuz;

         // Assign a random kinetic energy to the prompt neutron
         double ee = m_randomPromptEnergy();

         // Check if the neutron is outside of the bulk
         if( m_isOutside(xx, yy, zz) ) {
             int index = m_enter_other_section(section, xx, yy, zz);
             if (index >= 0) {  // >=0，因为组件数目在计算机中是从0开始计数
                 px[index][p_enter[index]] = xx;
                 py[index][p_enter[index]] = yy;
                 pz[index][p_enter[index]] = zz;
                 pe[index][p_enter[index]] = ee;
                 pux[index][p_enter[index]] = uux;
                 puy[index][p_enter[index]] = uuy;
                 puz[index][p_enter[index]] = uuz;
                 p_enter[index]++;
                 if (p_enter[index] >= m_nMax) {
                     p_enter[index] %= m_nMax;
                 }
                 number_of_enter[index]++;
             }
         } else{
                 px[section][p_enter[section]] = xx;
                 py[section][p_enter[section]] = yy;
                 pz[section][p_enter[section]] = zz;
                 pe[section][p_enter[section]] = ee;
                 pux[section][p_enter[section]] = uux;
                 puy[section][p_enter[section]] = uuy;
                 puz[section][p_enter[section]] = uuz;
                 p_enter[section]++;
                 if (p_enter[section] >= m_nMax) {
                     p_enter[section] %= m_nMax;
                 }
                 ng++;
             }
         } // end of prompt neutron loop
      } // end of isFission
            else if(isScattered) {
                // Get a random free path length for this energy
                  double d = m_randomPathLength(section);

                  // Send the scattered neutron to a isotropiaclly random direction in space
                  double phi  = 2.0*m_PI*m_randomFlat();
                  double theta= acos(2.0*m_randomFlat()-1.0);
                  double uux = sin(theta)*cos(phi);
                  double uuy = sin(theta)*sin(phi);
                  double uuz = cos(theta);
                  double xx  = x + d*uux;
                  double yy  = y + d*uuy;
                  double zz  = z + d*uuz;

                  // Scattering angle
                  // i.e. opening angle between new and old direction (it is found from dot product)
                  double cosThetaS = ux*uux + uy*uuy + uz*uuz;

                  // Determine the random kinetic energy to this scattered neutron
                  double ee = m_randomScatteredEnergy(e, cosThetaS, isElastic);
                  if(ee<0.0) isPossible = false;

                  // Check if the neutron is outside of the bulk
                  if( m_isOutside(xx, yy, zz) ) {
                      int index = m_enter_other_section(section, xx, yy, zz);
                      if (index >= 0) {
                          px[index][p_enter[index]] = xx;
                          py[index][p_enter[index]] = yy;
                          pz[index][p_enter[index]] = zz;
                          pe[index][p_enter[index]] = ee;
                          pux[index][p_enter[index]] = uux;
                          puy[index][p_enter[index]] = uuy;
                          puz[index][p_enter[index]] = uuz;
                          p_enter[index]++;
                          if (p_enter[index] >= m_nMax) {
                              p_enter[index] %= m_nMax;
                          }
                          number_of_enter[index]++;
                      }
                  } else if(isPossible){
                      px[section][p_enter[section]] = xx;
                      py[section][p_enter[section]] = yy;
                      pz[section][p_enter[section]] = zz;
                      pe[section][p_enter[section]] = ee;
                      pux[section][p_enter[section]] = uux;
                      puy[section][p_enter[section]] = uuy;
                      puz[section][p_enter[section]] = uuz;
                      p_enter[section]++;
                      if (p_enter[section] >= m_nMax) {
                          p_enter[section] %= m_nMax;
                      }
                       ng++;
                  }
                  else {
                      (void)0;
                  }
               } else if(isCaptured)
               {
                     (void)0;
               }
            return ng;
}

void NeutronTransport3D::start(){
   cout << endl;
   cout << "NeutronTransport 2.0 - May 2015" << endl;
   cout << "A Simplified Monte Carlo Neutron Transport Simulation " << endl << endl;

   cout << "--- Neutron Transport will start with the following initial values -------" << endl;
   cout << "Number of events        : " << fixed << setw(8) << setprecision(3) << m_evtMax << endl;
   cout << "Radius of the bulk (cm) : " << fixed << setw(8) << setprecision(3) << m_radius << endl;
   cout << "Seed of random # gen.   : " << setw(8) << m_randomSeed << endl;
   cout << "--------------------------------------------------------------------------" << endl;

   //*** Get cross-section data from the data file ***
   m_getData();

   cout << endl << "Press [return] to continue..." << endl;
   cout << "Here we utilize automatic pass this human enter!" << endl;
   // getchar();
}

//------------------------------------------------------------------------------------------------------------
// Calculates the reaction probabilities using neutron cross-section data evaluated
// for U235 and U238 for the given kinetic energy.
//------------------------------------------------------------------------------------------------------------
void NeutronTransport3D::m_getProbabilities(double energy)
{
   //for(int i=0; i<5; i++) m_micro235[i] = m_micro238[i] = 0.0;

    for (int s=0; s < ZONE; s++) {
       //=== Microscopic cross-sections ===
       if( energy <= m_eThermal ){  // 热中子的截面数据采用统一近似值
          m_micro235[s][0]=99.0;    m_micro238[s][0]=0.0;    // 0 -> Radiation capture
          m_micro235[s][1]=582.0;   m_micro238[s][1]=0.0;    // 1 -> Fission
          m_micro235[s][2]=13.78;   m_micro238[s][2]=8.871;  // 2 -> Elastic scattering
          m_micro235[s][3]=0.2;     m_micro238[s][3]=0.0;    // 3 -> In-elastic scattering
          m_micro235[s][4]=695.0;   m_micro238[s][4]=11.551; // 4 -> Total cross-sections
       }
       else{
          // values are obtained by linear extrapolation,线性插值获得快中子的截面数据
          for(int i=0; i<19; i++){
             if(energy >= m_Energy[i] && energy <= m_Energy[i+1]){
                double de = (energy-m_Energy[i])/(m_Energy[i+1]-m_Energy[i]);
                for(int j=0; j<5; j++){
                   m_micro235[s][j] = (m_s235[j][i+1]-m_s235[j][i])*de + m_s235[j][i];
                   m_micro238[s][j] = (m_s238[j][i+1]-m_s238[j][i])*de + m_s238[j][i];
                }
             }
          }
       }
    }

    for (int s=0; s<ZONE; s++) {
      //=== Macroscopic cross-sections ===
      for(int i=0; i<5; i++){
          m_macro235[s][i] = m_N235[s] * m_micro235[s][i] * m_barn;
          m_macro238[s][i] = m_N238[s] * m_micro238[s][i] * m_barn;
          m_macroAll[s][i] = m_macro235[s][i] + m_macro238[s][i];
      }
    }

  //=== mean free path of a neutron in the bulk ===, 计算获得平均自由程
    for (int s=0; s<ZONE; s++) {
        m_lambda[s] = 1.0/m_macroAll[s][4];
    }

   for (int s=0; s<ZONE; s++) {
      //=== Reaction probabilities ===
      m_pAbs[s] = (m_macroAll[s][0]+m_macroAll[s][1])/m_macroAll[s][4]; // absorbtion    probability wrt total cs
      m_pCap[s] = m_macroAll[s][0]/(m_macroAll[s][0]+m_macroAll[s][1]); // capture       probability wrt absorbtion cs
      m_pFis[s] = m_macroAll[s][1]/(m_macroAll[s][0]+m_macroAll[s][1]); // fission       probability wrt absorbtion cs
      m_p235[s] = m_macro235[s][1]/(m_macro235[s][1]+m_macro238[s][1]); // U235          probability wrt fission cs

      m_pSct[s] = (m_macroAll[s][2]+m_macroAll[s][3])/m_macroAll[s][4]; // scattering    probability wrt total cs
      m_pEla[s] = m_macroAll[s][2]/(m_macroAll[s][2]+m_macroAll[s][3]); // elastic s.    probability wrt scattering cs
      m_pInE[s] = m_macroAll[s][2]/(m_macroAll[s][2]+m_macroAll[s][3]); // in-elastic s. probability wrt scattering cs
   }
}

//------------------------------------------------------------------------------------------------------------
// Reads data from the 'fission-cross-section.data' to get energies (eV) and corresponding
// microscopic cross-sections (in barn) for U235 and U238.
// Datum are stroed in the arrays of m_Energy[j], m_s235[i][j] and m_s238[i][j] (i=0,4; j=0,19)
//
// ** You can download data file from:
//    http://www1.gantep.edu.tr/~bingul/simulation/fission/fission-cross-section.data
//
// ** Data is taken from the following web sites:
//
//    - http://www.ncnr.nist.gov/resources/n-lengths/list.html
//      Nuclear Center for Nuclear Resarch
//
//    - http://www.nndc.bnl.gov/
//      National Nuclear Data Center (Brookhaven National Lab.)
//      ENDF Evaluated Nuclear (reaction) Datas File is used.
//------------------------------------------------------------------------------------------------------------
void NeutronTransport3D::m_getData(void)
{
  ifstream csDataFile ("cross-section.data");

  cout << "Getting cross-section data... " << endl;

  if( csDataFile.is_open() )
  {
        // get 20-Kinetic energy value (in eV) for the cross-section calculation
        for(int j=0; j<20; j++)
          csDataFile >> m_Energy[j];

        // get U235 cross section (in barn) data
        for(int j=0; j<20; j++){
        for(int i=0; i<5;  i++)
          csDataFile >> m_s235[i][j];
        }

        // get U238 cross section (in barn) data
        for(int j=0; j<20; j++){
        for(int i=0; i<5;  i++)
          csDataFile >> m_s238[i][j];
        }
        csDataFile.close();
        cout << "done." << endl;
  }
  else{
   cout << "Cannot open 'fission-cross-section.data'." << endl;
   cout << "You can download it from:" << endl;
   cout << "http://www1.gantep.edu.tr/~bingul/simulation/fission/fission-cross-section.data" << endl;
   exit(1);
  }

  if( m_showDetails == false ) return;

  cout << "--- U235 cross-section data ----------------------------------------------" << endl;
  cout<< "  En. (ev)   Capture   Fission   Elestic   In-Elas   Total  " << endl;
  cout<< "  ========   =======   =======   =======   =======   =======" << endl;
  for(int j=0; j<20; j++){
     cout << scientific << setw(10) << setprecision(2) << m_Energy[j];
     for(int i=0; i<5;  i++)
        cout << fixed << setw(10) << setprecision(5) << m_s235[i][j];
  cout << endl;
  }
  cout << "--- U238 cross-section data ----------------------------------------------" << endl;
  cout<< "  En. (ev)   Capture   Fission   Elestic   In-Elas   Total  " << endl;
  cout<< "  ========   =======   =======   =======   =======   =======" << endl;
  for(int j=0; j<20; j++){
     cout << scientific << setw(10) << setprecision(2) << m_Energy[j];
     for(int i=0; i<5;  i++)
        cout << fixed << setw(10) << setprecision(5) << m_s238[i][j];
  cout << endl;
  }
  cout << "--------------------------------------------------------------------------" << endl << endl;
}

bool NeutronTransport3D::m_isOutside(double xx, double yy, double zz)
{
    if( (xx>-10 && xx<10.0) && (yy>-10 && yy<10.0) && (zz>-50 && zz<50.0) ) return false;
    return true;
}

//------------------------------------------------------------------------------------------------------------
// Returns nearest integer of a the double x
//------------------------------------------------------------------------------------------------------------
int NeutronTransport3D::m_nint(double x)
{
  if( (x - int(x)) >= 0.5) return int(x+1.0);
  else                     return int(x);
}


//***********************************************************************************************************
// Following methods, that performs MC simulation, use random numbers.
//
// m_randomFlat()            : Returns a uniform random deviate between 0.0 and 1.0.
// m_randomNeutronNumber()   : Returns a random neutron number after fission
// m_randomPathLength()      : Returns a random free path for a neutron (avr = m_lambda)
// m_randomPromptEnergy()    : Returns a random energy for a prompt neutron. (avr ~ 2.5 MeV)
// m_randomScatteredEnergy() : Returns a random energy for a scattered neutron.
//

//------------------------------------------------------------------------------------------------------------
// Returns a uniform random deviate between 0.0 and 1.0.
// Based on: Park and Miller's "Minimal Standard" random number generator (Comm. ACM, 31, 1192, 1988)
//------------------------------------------------------------------------------------------------------------
double NeutronTransport3D::m_randomFlat()
{
   const  int   im = 2147483647, ia = 16807;
   const  int   iq = 127773,     ir = 2836;
   const  double m = 128.0/im;
   int    k;
   double r;

   k = m_randomSeed / iq;
   m_randomSeed = ia*(m_randomSeed-k*iq) - ir*k;

   if(m_randomSeed < 0) m_randomSeed += im;

   r = m * (m_randomSeed/128);

   return r;
}

//------------------------------------------------------------------------------------------------------------
// Returns a normally distributed integer number of neutron generated per fission.
// The width of the distribution is sigma and mean is energy dependent.
// The result Nubar dpends also the nuclei.
//------------------------------------------------------------------------------------------------------------
int NeutronTransport3D::m_randomNubar(double energy, bool isU235)
{
  double mean, sigma = 1.0;
  double x, y, z;

  while(1){
        x = 2.0 * m_randomFlat() - 1.0;
        y = 2.0 * m_randomFlat() - 1.0;
        z = x*x + y*y;
        if( z <= 1.0 ) break;
  }

  // convert energy in eV to MeV
  energy *= 1.0e-6;

  // U235 nucleus
  if(isU235){
     if(energy<5.0) mean = 0.12*energy + 2.4;
     else           mean = 0.16*energy + 2.2;
  }
  // U238 nucleus
  else{
     if(energy<3.0) mean = 0.0666*energy + 2.4;
     else           mean = 0.1654*energy + 2.1538;
  }

  return m_nint(mean + sigma*x*sqrt(-2.0*log(z)/z));
}

//------------------------------------------------------------------------------------------------------------
// Returns a random path length obeying exponential decay law.
// This is actually a poission distribution exp(-x/lambda).
//------------------------------------------------------------------------------------------------------------

double NeutronTransport3D::m_randomPathLength(int section)
{
   return -m_lambda[section] * log(m_randomFlat());
}

//------------------------------------------------------------------------------------------------------------
// Returns a random kinetic energy in eV for a prompt neutron.
// The values are taken from a Maxwellian (or Watt) distribution.
// This distribution is still commonly used to describe the prompt fission neutron spectra.
//------------------------------------------------------------------------------------------------------------
double NeutronTransport3D::m_randomPromptEnergy()
{
  double emin =  0.05;      // minimum value of the kinetic energy in MeV
  double emax = 20.00;      // maximum value of the kinetic energy in MeV
  double Tm   =  1.29;      // Maxwellian temperature, somewhere between [1.290, 1.426]
  double eopt =  Tm/2.0;    // optimum value of the energy
  double scal = 2.0/(sqrt(m_PI)*pow(Tm,1.5));  // a scale factor
  double fmax = scal*sqrt(eopt)*exp(-eopt/Tm); // optimum value of the distribution

  // A rejection algorithm
  while(1)
  {
    double ptest = m_randomFlat()*fmax;
    double e     = m_randomFlat()*(emax-emin) + emin;
    double fmb   = scal*sqrt(e)*exp(-e/Tm);
    if (fmb > ptest) return e*1.0e+6;
  }
}

//------------------------------------------------------------------------------------------------------------
// Returns a kinetic energy in eV of the scattered neutron for the given
//           Ei : kinetic energy of the incident neutron in eV
// cosThetaScat : the cosine of scattering angle
//      elastic : the scattering type (true for elastic and false for in-elastic scattering).
// Note that aproximated values are used for the calculation.
//------------------------------------------------------------------------------------------------------------
double NeutronTransport3D::m_randomScatteredEnergy(double Ei, double cosThetaScat, bool elastic)
{
    double A = 237.0, Ef, Q;

    if(elastic){
        Ef = Ei*(A+cosThetaScat)*(A+cosThetaScat) / ((A+1.)*(A+1.));
        //Ef = Ei*(0.99 - m_randomFlat()/100.);
    }
    else{
        // Reaction Q value of the in-elastic reaction
        if(m_randomFlat()<0.1) Q = 148.0e+3;
        else                   Q =  44.9e+3;

        if(Q>Ei) Ef = -1.0; // reaction is not possible
        else     Ef = pow(sqrt(Ei)*cosThetaScat-A*sqrt(Ei-Q),2.) / ((A+1.)*(A+1.));
        //else     Ef = Ei*(0.99 - m_randomFlat()/100.) - Q;
        if(Ef<5.0e+04) Ef = -1.0; // no cross-section info
    }

    return Ef;
}

#endif /* FISSION_H_ */
