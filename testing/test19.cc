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

#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include "Fission.h"

#define PI  3.1415926

#define TRACE_PRINT(fmt, args...) do { \
                fprintf(stderr, "DEBUG: %s:%d:%s(): " fmt, \
                __FILE__, __LINE__, __func__, ##args); \
            } while(0)

using namespace std;

double Now()
{
    struct timespec timer;
    double tm;

    clock_gettime(CLOCK_REALTIME, &timer);
    tm = timer.tv_sec + timer.tv_nsec*0.000000001;

    return tm;
}

// Pass n lines
void File_FilterLine(ifstream &infile, int n)
{
    int i;
    string line;

    if (infile == NULL) {
        return;
    }

    for (i=0; i<n; i++) {
        if (infile.eof()) {
            break;
        }
        getline(infile, line);
#ifdef _DEBUG
        cout<<"File filter: "<<line<<endl;
#endif
    }
}

// string to double
// ENDF data lack of scientific mark e or E
double S2double(string s)
{
    int i, k;
    double a = 0.0;
    char buf[128] = {'\0'};
    char str[128] = {'\0'};

    strcpy(str, s.c_str());

    for (k=0, i=0; i<128; i++, k++) {
        if (str[k] == '\0') {
            buf[i] = '\0';
            break;
        }
        if (str[k] == '+' || str[k] == '-') {
            buf[i++] = 'e';
        }
        buf[i] = str[k];
    }

    a = strtod(buf, NULL);

    return a;
}

// Basic neutron cross-section class, read data from ENDF file
// type illustration:
//  31£¬(n, total)           MF=3 MT=1
//  32£¬(n, elastic)         MF=3 MT=2
//  318£¬(n, fission)        MF=3 MT=18
//  3102£¬(n, r)             MF=3 MT=102
//  1452, MF=1 MT=452       number of neutrons during every fission
//
class Nuclear_CrossSection {
private:
    int type;       // ENDF: type=31 represent MF=3 and MT=1, the total cross-section
    int np;         // number of cross-section datas
    double *energy;
    double *cross_section;
    void Print_CS3();
public:
    Nuclear_CrossSection(int type, std::string file);
    ~Nuclear_CrossSection();
    void Print_CrossSection();
};

const int CS_NUM_OF_EVERYLINE = 3;  // every line has how many effective datas
const int CS_SIGNIFICANT = 8;       // ENDF index width of scientific number

// file is the name of CS datas
Nuclear_CrossSection::Nuclear_CrossSection(int type, string file)
{
    ifstream infile;
    string line;
    char buffer[128] = {0};
    int i, k, n, tmp, line_of_data;

    this->type = type;
    energy = NULL;
    cross_section = NULL;
    infile.open(file.c_str(), ios::in);
    // check the file readable
    if (infile == NULL || infile.peek() == EOF) {
        cout<<"file "<<file<<" is not exist or empty!"<<endl;
        infile.close();
        return;
    }
    if (type == 31 || type == 32 || type == 318 || type == 3102) {
        getline(infile, line);
        sscanf(line.c_str(), "%d %d", &n, &tmp);
        np = n;

        if (n <= 0) {
            cout<<"file "<<file<<" store no cross section data!"<<endl;
            return;
        }
        energy = new double[n];
        cross_section = new double[n];

        line_of_data = (n/CS_NUM_OF_EVERYLINE) + 1;  // every line has 3 datas
        File_FilterLine(infile, 3);

        stringstream ss;
        string pa, pb;

        for ( n=0, i=1; i<line_of_data; i++) {
            getline(infile, line);
            ss << line;
            for (k=1; k<=CS_NUM_OF_EVERYLINE; k++) {
                ss >> pa >> pb;  // (E, cross-section)
                energy[n] = S2double(pa);
                cross_section[n] = S2double(pb);
                n++;
            }
            ss.clear();
        }
        // last line is special
        getline(infile, line);
        ss << line;
        tmp = np-(line_of_data-1)*CS_NUM_OF_EVERYLINE;
#ifdef _DEBUG
        cout<<"total lines are "<<line_of_data<<" , pairs of last line = "<<tmp<<endl;
#endif
        for (k=1; k<=tmp; k++) {
            ss>>pa>>pb;
            energy[n] = S2double(pa);
            cross_section[n] = S2double(pb);
#ifdef _DEBUG
            cout<<"last line string read : "<<pa<<" "<<pb<<endl;
            printf("double format : %.8le %.8le\n", energy[n], cross_section[n]);
#endif
            n++;
        }
        ss.clear();
#ifdef _DEBUG
        cout<<"Last line ("<<file<<"): "<<line<<endl;
        cout<<"NP = "<<np<<endl;
#endif
    }

    infile.close();
}

Nuclear_CrossSection::~Nuclear_CrossSection()
{
    if (cross_section != NULL) {
        delete[] cross_section;
    }
    if (energy != NULL) {
        delete[] energy;
    }
}

// Print MF=3 and MT=1 cross section data
void Nuclear_CrossSection::Print_CS3()
{
    int i;

    printf("  ev  \t  barns\n");
    printf("-------------------------------\n");
    for (i=0; i<np; i++) {
        printf("%.8le\t%.8le\n", energy[i], cross_section[i]);
    }
}

void Nuclear_CrossSection::Print_CrossSection()
{
    if (type == 31) {
        cout<<"MF=3, MT=1, total microscopic cross section :"<<endl;
        cout<<np<<endl;
        Print_CS3();
    } else if (type == 32) {
        cout<<"MF=3, MT=2, microscopic eslatic cross section :"<<endl;
        cout<<np<<endl;
        Print_CS3();
    } else if (type == 318) {
        cout<<"MF=3, MT=18, microscopic (n,f) cross section :"<<endl;
        cout<<np<<endl;
        Print_CS3();
    } else if (type == 3102) {
        cout<<"MF=3, MT=102, microscopic (n,r) cross section :"<<endl;
        cout<<np<<endl;
        Print_CS3();
    } else {
        (void)0;
    }
}

void Reaction_info()
{
    cout<<"1 barn = "<<1e-24<<" cm*cm"<<endl;
    cout<<"Microscopic cross section : [cm2], [barns] or [m2] unit"<<endl;
    cout<<"Macroscopic cross section : [cm-1] or [m-1] unit"<<endl;
}

void MonteCarlo_nt()
{
    double start, end;
   int   maxEvent = 10000;  // Maximum number of events (number of neutrons) required (default  1)
   int     maxGen =   1000;  // Maximum number of generations (defalult 100)
   double  purity = 0.01;  // Purity of 235 in the medium (default  3%)
   double  radius;         // Radius of the spherical bulk (default  10 cm)
   bool    detail = false;  // To print details of each event and generations (default false)
   int     seed   = 31415;  // A seed to initiate random number generator (default 0)

   radius = 21.0;
   NeutronTransport nt(maxEvent, maxGen, purity, radius, detail, seed);

   start = Now();
   nt.start();
   nt.execute();
   nt.end();
   end = Now();
   cout<<"Total run time is "<<end-start<<" sec"<<endl;
}

int main(int argc, char *argv[])
{
    TRACE_PRINT("start to integrated example 1!\n");
    TRACE_PRINT("FOUR ASSEMBLIES : Monte Carlo NEUTRON SIMULATION!\n");

    Reaction_info();

    // read cross section parameters and print it
    Nuclear_CrossSection U235(3102, "U235_MF3_MT102.txt");
    U235.Print_CrossSection();

    MonteCarlo_nt();

    return 0;
}

