#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <vector>
#include <cmath> // For sqrt and pow
#include <set>
#include <map>
#include <iostream>
using namespace std;

// --- GEOMETRIC & ENTITY STRUCTURES ---

struct Point {
    double x, y;
};

// --- UTILITY FUNCTIONS ---

/**
 * @brief Calculates the Euclidean distance between two points.
 */
inline double distance(const Point& p1, const Point& p2) {
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}


// --- PROBLEM & SOLUTION STRUCTURES (remaining definitions are the same) ---

/**
 * @brief containa all the information related to the package
 * 
 * @param weight weight of the package
 * @param value value of the package
 * @param density value / weight
 * @param id id number of the package 0-dry 1-per 2-othe
 */
struct PackageInfo {
    double weight, value, density;
    int id;

    bool operator<(const PackageInfo& other) const {
        return density > other.density;
    }
};

struct Packages {
    int dry_food;
    int perishable_food;
    int other_supplies;
};

struct Village {
    int id;
    Point coords;
    int population;
    int rem_population_food;
    Packages pack_received;
    int rem_population_other;  // added people remaining to be served

    /**
     * @brief call this function at any instant and it will tell the value added to this 
     * 
     * @param packs give out the packet information here, make sure ordering is maintained
     * @return the value added to this village only when pack_received is filled
     * @note gives priority to perishables , then dry then others
     * 
     * @note make sure pack_received is populated correctly before using this fucntion
     */
    double value_receieved(vector<PackageInfo> packs) {
        double value = 0;

        int food_needed = population * 9;
        value += min(pack_received.perishable_food, food_needed) * packs[1].value;

        if (pack_received.perishable_food >= food_needed) return value;
        int rem =  food_needed - pack_received.perishable_food;
        value += min(rem, pack_received.dry_food) * packs[0].value;

        if (pack_received.dry_food >= 9*rem) return value;

        value += min(population, pack_received.other_supplies) * packs[2].value;

        return value;
        
    }

};
struct Drop {
    Village v;
    int village_id;
    int dry_food;
    int perishable_food;
    int other_supplies;
};

struct Trip {
    int dry_food_pickup;
    int perishable_food_pickup;
    int other_supplies_pickup;
    vector<Drop> drops;
};
struct Helicopter {
    int id;
    int home_city_id;
    double d_max;
    Point home_city_coords;
    double weight_capacity;
    double distance_capacity;
    double fixed_cost; // F
    double alpha;
    vector<Trip> trips;

    /**
     * @brief calculates the distance in a trip
     * 
     * @param t the trip for which the distance is to be calculated
     * 
     * @return distance travelled in a trip only
     */
    double trip_dist(Trip t) {
        
        double dist = 0;
        Point prev_coords = home_city_coords;
        for (Drop d: t.drops) {
            dist += distance(prev_coords, d.v.coords);
            prev_coords = d.v.coords;
        }
        dist += distance(prev_coords, home_city_coords);
        return dist;
    }

    /**
     * @brief calculates the distance across a vector of Trips
     * 
     * @param trips a vector of trips
     * 
     * @return sum of distance across each trips
     */
    double trips_dist(vector<Trip> _trips){

        double _trips_dist = 0;
        for(Trip t : _trips){
            _trips_dist += trip_dist(t);
        }
        return _trips_dist;
    }

    /**
     * @brief use thus for per trip cost
     * 
     * @param t trip for which the cost to find
     * 
     * @return one trip cost
     */
    double trip_cost(Trip t) {
        double _trip_dist = trip_dist(t);
        double _trip_cost = fixed_cost + alpha * _trip_dist;
        return _trip_cost;
    }

    /**
     * @brief finds the cost of trip across several cities
     * 
     * @param _trips vector of all the trips thaat are to be considered for the trip cost
     * 
     * @return total trip cost across all the trips given
     */
    double trips_cost(vector<Trip> _trips) {
        double cost = 0;
        for (Trip t: _trips) {
            cost += trip_cost(t);
        }
        return cost;
    }

    /**
     * @brief checks if the trip is within the specified limit
     * 
     * @param t trip to check the limit on
     * 
     * @return true if travel was smaller than dcap else false
     */
    bool check_dcap(Trip t) {
        double _trip_dist = trip_dist(t);
        // cout << _trip_dist << "&" << distance_capacity << endl;
        return distance_capacity >= _trip_dist;
    }

    /**
     * @brief check if all the trips obey the Dmax for the problem which is per Helicopter
     * 
     * @param trips all the vector of trips upon which the d_max cap is to be applied
     * 
     * @return true if d_max cap is satisfied else false
     */
    bool check_dmax(vector<Trip> _trips) {
        double total_trips_dist = trips_dist(_trips);
        return total_trips_dist <= d_max;
    }

    /**
     * @brief prepares a drop variable for return 
     * 
     * @return a Drop
     * @note no changes to v being done
     */
    Drop prepare_drop(Village v) {
        Drop d;
        d.dry_food = 0;
        d.other_supplies = 0;
        d.perishable_food = 0;
        d.v  = v;
        d.village_id = d.v.id;
        return d;
    }

    /**
     * @brief prepares a trip variable for return 
     * 
     * @return a Trip
     * @note no change done to v
     */
    Trip prepare_trip(Village v) {
        Trip t;
        Drop d = prepare_drop(v);
        t.drops = {d};
        t.dry_food_pickup = 0;
        t.other_supplies_pickup = 0;
        t.perishable_food_pickup = 0;

        return t;
    }

    vector<Trip> try_new_trip(map<int, Village> vmap) {
        vector<Trip> sugg;
        for (pair<int, Village> p : vmap) {
            vector<Trip> trips_for_h = trips;
            Trip t = prepare_trip(p.second);
            trips_for_h.push_back(t);

            if (check_dcap(t) && check_dmax(trips_for_h)) {
                sugg.push_back(t);
            }
        }

        return sugg;
    }

    /**
     * @brief check by adding only one village at a time and ensures it is valud for both dcap and dmax
     * 
     * @param t_idx index in the trips for ease of d_max
     * 
     * @return vector trips you can get out of only one single trip, then you can use these to wrap them with the orihinal set of of trips but this chamges and make that a child
     */
    vector<Trip> try_adding_village(int t_idx, map<int, Village> vmap) { 
        Trip t = trips[t_idx];
        vector<Trip> sug_trips;
        vector<Trip> trip_for_h = trips;
        // cout << "trips for h" << endl;
        // for (size_t m_idx =0; m_idx<trip_for_h.size(); m_idx++) {
        //     for (size_t v_idx=0; v_idx<trip_for_h[m_idx].drops.size(); v_idx++) {
        //         cout << trip_for_h[m_idx].drops[v_idx].village_id << " ";
        //     }
        //     cout << endl;
        // }
        // cout << "trips for h over" << endl;
        set<int> vids;
        // cout << "start" << endl;
        for (Drop d: t.drops) {
            // cout << "try adding "<< d.village_id << " " << d.v.id << endl;
            vids.insert(d.village_id);
        }
        // cout << endl;
        for (pair<int, Village> p: vmap) {
            if (vids.find(p.first)!=vids.end()) continue; // for uniqueness of the villages in a trip
            Village new_vill = p.second;
            cout << "new vill " << new_vill.id << endl;
            Trip sug_trip(t);
            Drop new_vil_drop = prepare_drop(new_vill);
            sug_trip.drops.push_back(new_vil_drop);

            trip_for_h[t_idx] = sug_trip;
            // perform check
            if (check_dcap(sug_trip) && check_dmax(trip_for_h)) {
                sug_trips.push_back(sug_trip);
            }
        }

        return sug_trips;
    }
};

struct ProblemData {
    double time_limit_minutes;
    double d_max;
    vector<PackageInfo> packages;
    vector<Point> cities;
    vector<Village> villages;
    vector<Helicopter> helicopters;
};




struct HelicopterPlan {
    int helicopter_id;
    Helicopter heli;

};

using Solution = vector<HelicopterPlan>;

#endif // STRUCTURES_H