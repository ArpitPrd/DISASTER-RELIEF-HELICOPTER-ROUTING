#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <vector>
#include <cmath> // For sqrt and pow
#include <unordered_map>
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

struct PackageInfo {
    double weight, value, density;
    int id;
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
     */
    double value_receieved(vector<PackageInfo> packs) {
        double value = 0;

        value += min(pack_received.perishable_food, population) * packs[1].value;

        if (pack_received.perishable_food >= population) return value;
        int rem1 =  population - pack_received.perishable_food;
        value += min(rem1, pack_received.dry_food) * packs[0].value;

        if (pack_received.dry_food >= rem1) return value;

        int rem2 = rem1-pack_received.dry_food;

        value += min(rem2, pack_received.other_supplies) * packs[2].value;

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


    bool check_dcap(double dcap, Trip t, Point home) {
        double trip_dist = trip_distance_travelled(t, home);
        return dcap >= trip_dist;
    }

    bool check_dmax(double dmax, vector<Trip> trips, Point home) {
        double total_trips_dist = all_trip_distance(trips, home);
        return total_trips_dist <= dmax;
    }

    vector<Trip> try_adding_village(Trip t) { 
        vector<Trip> sug_trips;
        for (pair<int, Village*> p: vmap) {

            Village new_vill = *p.second;

            Trip sug_trip(t);
            Drop new_vil_drop = prepare_drop(new_vill);
            sug_trip.drops.push_back(new_vil_drop);

            // perform check
            if (check_dcap(distance_capacity, sug_trip, home_city_coords) && check_dmax(d_max, trips, home_city_coords)) {
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

    /**
     * @brief calculates the distance in a trip
     * 
     * @param t the trip for which the distance is to be calculated
     * @param home the home coordinates for the center
     */
    double trip_dist(Trip t, Point home) {
        
        double dist = 0;
        Point prev_coords = home;
        for (Drop d: t.drops) {
            dist += distance(prev_coords, d.v.coords);
            prev_coords = d.v.coords;
        }
        dist += distance(prev_coords, home);
        return dist;
    }

    /**
     * @brief calculates the distance across a vector of Trips
     * 
     * @param trips a vector of trips
     * @param home point of home
     * 
     * @return sum of distance across each trips
     */
    float trips_dist(){
        Point home = heli.home_city_coords;

        float _trips_dist = 0;
        for(Trip t : heli.trips){
            _trips_dist += trip_dist(t, home);
        }
        return _trips_dist;
    }

    /**
     * @brief use thus for per trip cost
     * 
     * @param t trip for which the cost to find
     * @param home about which the trip happened
     * @param F fixed cost for the trip
     * @param alpha variable cost for the trip
     * 
     * @return one trip cost
     */
    double trip_cost(Trip t, Point home, double F, double alpha) {
        double _trip_dist = trip_dist(t, home);
        double _trip_cost = F + alpha * _trip_dist;
        return _trip_cost;
    }

    /**
     * @brief finds the cost of trip across several cities
     * 
     * @return total trip cost across all the trips given
     */
    double trips_cost() {
        double cost = 0;
        for (Trip t: heli.trips) {
            cost += trip_cost(t, heli.home_city_coords, heli.fixed_cost, heli.alpha);
        }
        return cost;
    }


    /**
     * @brief find the value added in exactly one trip
     */

    double trip_value(Trip t, vector<PackageInfo> packs) {
        double value = 0;
        for (Drop d: t.drops) {
            value += d.v.value_receieved(packs);
        }
    }

    /**
     * @brief returns the value added by all the trips
     */
    double trips_value(vector<PackageInfo> packs) {
        double value = 0;

        for (Trip t: heli.trips) {
            value += trip_value(t, packs);
        }
    }

};

using Solution = vector<HelicopterPlan>;

#endif // STRUCTURES_H