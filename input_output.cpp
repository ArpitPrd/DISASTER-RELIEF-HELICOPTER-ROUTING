#include "stdc++.h"

struct village{
    int xv, yv, nv;
};

struct city{
    int xc, yc;
};

struct packet{
    int type;
    float weight, value;
    // type = 0 -> wet
    // type = 1 -> dry
    // type = 2 -> other
};

struct helicopter{
    int city_id, wcap, dcap, F, alpha;
    // wcap -> weight capacity per trip
    // dcap -> distance capacity per trip
    // Dmax -> total distance capacity
};

int trip_cost(int dist, int F, int alpha){
    return F + alpha*dist;
}

int main(){
    // ----------------------------INPUT------------------------------------------------------------------
    int processing_time;  // in minutes
    cin >> processing_time;
    int Dmax;
    cin >> Dmax;
    packet wet, dry, other;
    wet.type = 0;
    dry.type = 1;
    other.type = 2;
    cin >> wet.weight >> wet.value >> dry.weight >> dry.value >> other.weight >> other.value;
    int n_cities;
    cin >> n_cities;
    vector<city> C(n_cities);
    for(int i = 0; i < n_cities; i++){
        cin >> C[i].xc >> C[i].yc;
    }
    int n_vilage;
    cin >> n_vilage;
    vector<village> V(n_vilage);
    for(int i= 0; i < n_vilage; i++){
        cin >> V[i].xv >> V[i].yv >> V[i].nv;
    }
    int n_h;
    cin >> n_h;
    vector<helicopter> H(n_h);
    for(int i = 0; i < n_h; i++){
        cin >> H[i].city_id >> H[i].wcap >> H[i].dcap >> H[i].F >> H[i].alpha;
    }
    //----------------------------------------------------------------------------------------------------

}
