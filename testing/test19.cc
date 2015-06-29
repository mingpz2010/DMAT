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
#include <cstdio>
#include <cstdlib>
#include <string>

#define PI  3.1415926

#define TRACE_PRINT(fmt, args...) do { \
                fprintf(stderr, "DEBUG: %s:%d:%s(): " fmt, \
                __FILE__, __LINE__, __func__, ##args); \
            } while(0)

using namespace std;

// ���������ӽ����࣬���ڴ��ļ���ȡ�������ݲ����д����ṩ������õ�����
// type (���ڱ�ʾ���ϵĽ�������) :
//  31��(n, total), MF=3 MT=1, �����ܵ�΢�۽���
//  32��(n, elastic), MF=3 MT=2�����ϵ�΢�۵���ɢ�����
//  318��(n, fission)��MF=3 MT=18�����ϵ�΢���ѱ����
//  3102��(n, r), MF=3 MT=102�����ϵ�΢��r����������
//  1452, MF=1 MT=452������number of neutrons during every fission
//
//
class Nuclear_CrossSection {
private:
    int type;  // ��������ͣ�����ENDF�Ķ��壬����type=31����MF=3,MT=1,��������ܵ�΢�۽���
    int np;  // �������ݵĸ���
    double *energy;
    double *cross_section;
    void Print_CS3();
public:
    Nuclear_CrossSection(int type, std::string file);
    ~Nuclear_CrossSection();
    void Print_CrossSection();
};

const int CS_NUM_OF_EVERYLINE = 3;  // ÿ�н��������м�����Ч����
const int CS_SIGNIFICANT = 8;       // ENDF��ָ����ʽ��︡�����ĵ������ֵ�ʵ����Чλ��

// ����fileΪ�����ļ�������
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
    // ��׳���ж�
    if (infile == NULL || infile.peek() == EOF) {
        cout<<"file "<<file<<" is not exist or empty!"<<endl;
        infile.close();
        return;
    }
    if (type == 31 || type == 32 || type == 318 || type == 3102) { // ΢�۽����ļ��İ��ж�ȡ
        getline(infile, line);
        sscanf(line.c_str(), "%d %d", &n, &tmp);
        np = n;

        if (n <= 0) {
            cout<<"file "<<file<<" store no cross section data!"<<endl;
            return;
        }
        energy = new double[n];
        cross_section = new double[n];

        line_of_data = (n/CS_NUM_OF_EVERYLINE) + 1;  // ÿ����������
        File_FilterLine(infile, 3);

        stringstream ss;
        string pa, pb;
        // ǰ���line_of_data-1�ж�ȡ, n��������Ե���ŵ���
        for ( n=0, i=1; i<line_of_data; i++) {
            getline(infile, line);
            ss << line;
            for (k=1; k<=CS_NUM_OF_EVERYLINE; k++) {
                ss >> pa >> pb;  // ÿ�������ȡ(����-����)��
                energy[n] = S2double(pa);
                cross_section[n] = S2double(pb);
                n++;
            }
            ss.clear();
        }
        // ���һ�еĶ�ȡ��Ҫ���⴦��
        getline(infile, line);
        ss << line;
        tmp = np-(line_of_data-1)*CS_NUM_OF_EVERYLINE;  // ���һ���ж��ٽ����
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

// ��ӡMF=3 MT=1�Ľ�������
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

// ��infileָ����ļ����ӵ�ǰָ��λ�������˵�n��
void File_FilterLine(ifstream &infile, int n)
{
    int i;
    string line;

    if (infile == NULL) {
        return;
    }

    // �����е�ʱ������������������ļ�β����ֱ�ӵ��ļ�β��
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

// ����string�ַ���תdouble����
// ����ENDF���ݿ������ݳ�Ϊ1.003235+7��ȱ�ٿ�ѧ�������е�e��E
// ���Ա����������e���ַ����У�������strtodת��Ϊdouble����
// (ע���ú�������ʹ�õ�ǰ����������ַ���s������һ��ENDF�е�double����)
double S2double(string s)
{
    int i, k;
    double a = 0.0;
    char buf[128] = {'\0'};
    char str[128] = {'\0'};

    strcpy(str, s.c_str());

    // kָ��Ϊs�ַ�����iָ��Ϊ���e��E���ַ���
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

void Reaction_info()
{
    cout<<"1 barn = "<<1e-24<<" cm*cm"<<endl;
    cout<<"Microscopic cross section : [cm2], [barns] or [m2] unit"<<endl;
    cout<<"Macroscopic cross section : [cm-1] or [m-1] unit"<<endl;
}

int main(int argc, char *argv[])
{
    TRACE_PRINT("start to integrated example 1!\n");
    TRACE_PRINT("FOUR ASSEMBLIES : Monte Carlo NEUTRON SIMULATION!\n");

    Reaction_info();

    // read cross section parameters and print it
    Nuclear_CrossSection U235(3102, "U235_MF3_MT102.txt");
    U235.Print_CrossSection();

    return 0;
}

