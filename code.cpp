#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
#define PI 3.1415
#define n 100

int u_final[15000];
double t_final[15000];
double m_final[15000];
double h_final[15000];
double v_final[15000];

class Data {
public:
    double t;
    double m;
    double phi;
    double psi;
    Data(double t, double m, double phi, double psi) :t(t), m(m), phi(phi), psi(psi) {}
};

vector<Data> vec;
//
void task_Koshe()
{
    double t, dt = 0.01, t_max = 100, mmax = 0;
    double old_psi_1, old_psi_2, old_psi_3;
    double old_h, old_m, old_v, u;
    double new_h, new_m, new_v, new_psi_1, new_psi_2, new_psi_3;
    for (double phi = 0; phi <= 2 * PI; phi += 2 * PI / n)
    {
        for (double psi = -PI / 2; psi <= PI / 2; psi += PI / n)
        {
            t = 0;
            old_psi_1 = 1000 * cos(phi) * cos(psi), old_psi_2 = 1000 * sin(phi) * cos(psi), old_psi_3 = 1000 * sin(psi);
            old_h = 0, old_m = 20000, old_v = 0, u = 0;
            new_h = 0, new_m = 0, new_v = 0, new_psi_1 = 0, new_psi_2 = 0, new_psi_3 = 0;
            while (t <= t_max) {
                if (old_psi_3 * 3000 / old_m - old_psi_2 < 0 || old_m <= 2000) {
                    u = 0;
                }
                else {
                    u = 500;
                }
                new_m = old_m - u * dt;
                new_v = old_v + (-9.8 + u * 3000 / old_m - 1.25 * exp(-0.00013 * old_h) * old_v * old_v / old_m) * dt;
                new_h = old_h + (old_v)*dt;
                new_psi_1 = old_psi_1 + (-0.00013 * old_psi_3 * 1.25 * exp(-0.00013 * old_h) * old_v * old_v / old_m) * dt;
                new_psi_2 = old_psi_2 + old_psi_3 * ((u * 3000 / (old_m * old_m)) - (1.25 * exp(-0.00013 * old_h) * old_v * old_v / (old_m * old_m))) * dt;
                new_psi_3 = old_psi_3 + (((2 * old_psi_3 * 1.25 * exp(-0.00013 * old_h) * old_v) / old_m) - old_psi_1) * dt;
                t += dt;
                if (new_h >= 50000 && new_v >= 1500)
                {
                    vec.push_back(Data(t, new_m, phi, psi));
                   // if (new_m > mmax) mmax = new_m;
                    break;
                }
                if (new_v <= 0)
                    break;
                old_m = new_m ;
                old_v = new_v ;
                old_h = new_h;
                old_psi_1 = new_psi_1;
                old_psi_2 = new_psi_2 ;
                old_psi_3 = new_psi_3;
            }

        }

    }
    //cout << mmax << endl;
}
void make_grapic(double phi, double psi) 
{
    int count = 0;
    double t=0, dt = 0.01, t_max = 100,  mmax=0;
    double old_psi_1 = cos(phi) * cos(psi), old_psi_2 = sin(phi) * cos(psi),old_psi_3 = sin(psi);
    double old_h = 0, old_m = 20000, old_v = 0, u = 0;
    double new_h = 0, new_m = 0, new_v = 0, new_psi_1 = 0, new_psi_2 = 0, new_psi_3 = 0;
            while (t <= t_max) {
                if (old_psi_3 * 3000 / old_m - old_psi_2 < 0 || old_m <= 2000) {
                    u = 0;
                }
                else {
                    u = 500;
                }
                new_m = old_m - u * dt;
                new_v = old_v + (-9.8 + u * 3000 / old_m - 1.25 * exp(-0.00013 * old_h) * old_v * old_v / old_m) * dt;
                new_h = old_h + (old_v)*dt;
                new_psi_1 = old_psi_1 + (-0.00013 * old_psi_3 * 1.25 * exp(-0.00013 * old_h) * old_v * old_v / old_m) * dt;
                new_psi_2 = old_psi_2 + old_psi_3 * ((u * 3000 / (old_m * old_m)) - (1.25 * exp(-0.00013 * old_h) * old_v * old_v / (old_m * old_m))) * dt;
                new_psi_3 = old_psi_3 + (((2 * old_psi_3 * 1.25 * exp(-0.00013 * old_h) * old_v) / old_m) - old_psi_1) * dt;
                
                u_final[count] = u;
                t_final[count] = t;
                m_final[count] = old_m;
                h_final[count] = old_h;
                v_final[count] = old_v;

                t += dt;
                count++;
                old_m = new_m;
                old_v = new_v;
                old_h = new_h;
                old_psi_1 = new_psi_1;
                old_psi_2 = new_psi_2;
                old_psi_3 = new_psi_3;
            }

       
}


int main()
{
    task_Koshe();
    double m_max = vec[0].m;
    Data res = vec[0];
    for (Data d : vec)
    {
        if (d.m > m_max) {
            m_max = d.m;
            res = d;
        }
        //cout<<d.t<<" "<< d.m <<" " << d.phi << " " << d.psi << endl;
    }
    cout << endl;
    cout <<"t= " << res.t << endl << "m= " << res.m << endl << "phi= " << res.phi << endl << "psi= " << res.psi << endl << endl;
    make_grapic(res.phi,res.psi);
    for (int i = 0; i < 10000; i++) {
        if (i <= res.t / 0.01)
            cout <<"t= " << t_final[i] << "    " << "u= " << u_final[i] << "    " << "m= " << m_final[i] << "    " <<"h= " << h_final[i] << "    " <<"v= " << v_final[i] << endl;
    }
        return 0;
}
