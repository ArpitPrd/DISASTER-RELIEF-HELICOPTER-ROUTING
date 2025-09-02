#include "solver.h"
#include <iostream>
#include <chrono>
#include <map>
#include <algorithm>

using namespace std;
using namespace std::chrono;
using hrc = time_point<high_resolution_clock>;
#define now() high_resolution_clock::now()
// You can add any helper functions or classes you need here.

map<int, Helicopter*> hmap;
map<int, Village*> vmap;

/**
 * @brief calculates the distance in a trip
 */
float trip_distance_travelled(Trip t, Helicopter h) {
    float dist = 0;
    Point prev_coords = h.home_city_coords;
    for (Drop d: t.drops) {
        dist += distance(prev_coords, d.v.coords);
    }
    return dist;
}

/**
 * @brief use thus for per trip cost
 */
float trip_cost(Trip t, Helicopter h) {
    float dist = trip_distance_travelled(t, h);
    float cost = h.fixed_cost + h.alpha*dist;
    return cost;
}

/**
 * @brief use this for calculation of scores for one trip
 */
float all_trip_cost(HelicopterPlan hp) {
    Helicopter h = *hmap[hp.helicopter_id];
    float sum_across_trips = 0;
    for (Trip trip : h.trips) {
        sum_across_trips += trip_cost(trip, h);
    }
    return sum_across_trips;
}

float plan_cost(vector<HelicopterPlan> hps) {
    float cost = 0;
    for (auto hp: hps) {
        cost += all_trip_cost(hp);
    }
}

/**
 * @brief general purpose Timer contains term criteria, restart criteria
 * 
 * @param start_ts abs start time
 * @param end_time relative ending time
 * @param last_restart_time abs last restart time
 * @param eps_restart rel time lapse after which timer indicates the solver for a restart
 */
class Timer {
public:
    hrc start_ts;
    int end_time;
    hrc last_restart_ts;
    int eps_restart;
    Timer(int end_time, int eps_restart) {
        this->end_time = end_time;
        this->eps_restart = eps_restart;
    }
    void start() {
        start_ts = now();
    }
    bool check_term() {
        hrc ts = now();
        auto duration = duration_cast<seconds>(ts-start_ts);
        if (duration.count() > end_time + 1) {
            return true;
        }
        return false; 
    }
    float get_time() {
        hrc ts = now();
        auto duration = duration_cast<seconds>(ts-start_ts);
        return duration.count();
    }
    bool restart() {
        hrc ts = now();
        auto duration = duration_cast<seconds>(ts-last_restart_ts);
        last_restart_ts = ts;
        if (duration.count() > eps_restart) {
            return true;
        }
        return false;
    }
};

/**
 * @brief get all villages within distance dcap from center village
 * 
 * @param villages list of all villages
 * @param center village from which we want to find other villages within distance dcap
 * @param dcap distance capacity
 * @return std::vector<Village> list of villages within distance dcap from center village
 */
std::vector<Village> villagesWithinDistance(
    const std::vector<Village>& villages,
    double dcap,
    const Point& p
) {
    std::vector<Village> result;
    for (const auto& v : villages) {
        if (distance(v.coords, p) <= dcap) {
            result.push_back(v);
        }
    }
    return result;
}

/**
 * @brief get all possible paths (up to 2 villages) within distance d from point p
 * @note returns trips in sorted order of distance (shortest first)
 * 
 * @param p point from which we want to find paths
 * @param villages list of all villages
 * @param d maximum distance
 * @return std::vector<std::vector<Village>> list of all possible paths (up to 2 villages) within distance d from point p
 */
std::vector<Trip> planAllPaths(
    const Point& p, const std::vector<Village>& villages, double d) 
{
    std::vector<Trip> allPaths;

    // Single-village paths
    for (const auto& v : villages) {
        double totalDist = distance(p, v.coords) + distance(v.coords, p);
        if (totalDist <= d) {
            Drop drop;
            drop.village_id = v.id;
            drop.v = v;
            drop.dry_food = 0;
            drop.perishable_food = 0;
            drop.other_supplies = 0;
            Trip trip;
            trip.dry_food_pickup = 0;
            trip.perishable_food_pickup = 0;
            trip.other_supplies_pickup = 0;
            trip.drops.push_back(drop);
            trip.trip_distance = totalDist;
            allPaths.push_back(trip);
        }
    }

    // Two-village paths
    for (size_t i = 0; i < villages.size(); i++) {
        for (size_t j = i + 1; j < villages.size(); j++) {
            double totalDist = distance(p, villages[i].coords)
                             + distance(villages[i].coords, villages[j].coords)
                             + distance(villages[j].coords, p);
            if (totalDist <= d) {
                Drop drop1;
                drop1.village_id = villages[i].id;
                drop1.v = villages[i];
                drop1.dry_food = 0;
                drop1.perishable_food = 0;
                drop1.other_supplies = 0;   
                Drop drop2;
                drop2.village_id = villages[j].id;
                drop2.v = villages[j];
                drop2.dry_food = 0;
                drop2.perishable_food = 0;
                drop2.other_supplies = 0;
                Trip trip;
                trip.dry_food_pickup = 0;
                trip.perishable_food_pickup = 0;
                trip.other_supplies_pickup = 0;
                trip.drops.push_back(drop1);
                trip.drops.push_back(drop2);
                trip.trip_distance = totalDist;
                allPaths.push_back(trip);
            }
        }
    }

    // Sort by total distance
    std::sort(allPaths.begin(), allPaths.end(),
              [](const Trip& a, const Trip& b) {
                  return a.trip_distance < b.trip_distance;
              });


    return allPaths;
}

/**
 * @brief select trips (sorted in ascending order) from the list such that their total distance is within Dmax
 * 
 * @param trips list of all possible trips (sorted by distance)
 * @param Dmax maximum total distance
 * @return std::vector<Trip> list of selected trips within distance Dmax
 */
vector<Trip> selectTripsWithinLimit(const vector<Trip>& trips, double Dmax) {
    vector<Trip> result;
    double totalDist = 0.0;

    for (const auto& trip : trips) {
        if (totalDist + trip.trip_distance <= Dmax) {
            result.push_back(trip);
            totalDist += trip.trip_distance;
        } else {
            break; // stop as soon as the budget is exceeded
        }
    }
    return result;
}

/**
 * @brief assign packages to a trip uniformly across all villages in the trip
 * 
 * @param t trip to which packages are to be assigned
 * @param wcap weight capacity of the helicopter
 * @param pack list of package types with their weights and values
 * @return Trip updated trip with assigned packages
 */
Trip assignPackages(Trip t, int wcap, vector<PackageInfo>& pack) {
    int numVillages = (int)t.drops.size();
    int totalAssigned = 0;

    // priority order: perishable (1), dry (0), other (2)
    vector<int> order = {1, 0, 2};

    for (int pkgType : order) {
        int remainingCap = wcap - totalAssigned;
        if (remainingCap <= 0) break;

        // find maximum integer X satisfying constraint
        int X = remainingCap / (numVillages * pack[pkgType].weight);
        if (X <= 0) continue;

        int assignedThisType = 0;

        // Assign uniformly across all villages
        for (auto& d : t.drops) {
            if (pkgType == 0) d.dry_food = X;
            if (pkgType == 1) d.perishable_food = X;
            if (pkgType == 2) d.other_supplies = X;
            assignedThisType += X * pack[pkgType].weight;
        }

        // Adjust if overshoot
        while (assignedThisType > remainingCap) {
            for (auto& d : t.drops) {
                if (assignedThisType <= remainingCap) break;
                if (pkgType == 0 && d.dry_food > 0) {
                    d.dry_food -= 1;
                    assignedThisType -= pack[pkgType].weight;
                }
                if (pkgType == 1 && d.perishable_food > 0) {
                    d.perishable_food -= 1;
                    assignedThisType -= pack[pkgType].weight;
                }
                if (pkgType == 2 && d.other_supplies > 0) {
                    d.other_supplies -= 1;
                    assignedThisType -= pack[pkgType].weight;
                }
            }
        }

        totalAssigned += assignedThisType;
    }

    // update trip totals
    t.dry_food_pickup = 0;
    t.perishable_food_pickup = 0;
    t.other_supplies_pickup = 0;
    for (auto& d : t.drops) {
        t.dry_food_pickup += d.dry_food;
        t.perishable_food_pickup += d.perishable_food;
        t.other_supplies_pickup += d.other_supplies;
    }

    return t;
}




/**
 * @brief Defines the search problem space
 * 
 * @param b_states best states
 * @param data contains the information how each state is going to look like, this has been put in the Space because it contains information about how states are ging to look like and the sampler needs to know this
 * @def add_to_lm: add this "state" to local maxima
 * @def add_to_p: add this "state" to potential searches, can remove this later
 * 
 * @remark sampler is put here since you "sample" from a "space"
 */
class HCSpace {
public:
    vector<HCState> b_states;
    
    ProblemData data;
    bool red;
    
    HCSpace(ProblemData _data, bool _red) {
        data = _data;
        red = _red;
    }

    HCState sample() {
        vector<HelicopterPlan> hplans;
        for (const auto& heli : data.helicopters) {
            HelicopterPlan hplan;
            hplan.helicopter_id = heli.id;
            // get all villages within distance capacity from home city
            auto nearby_villages = villagesWithinDistance(data.villages, heli.distance_capacity, data.cities[heli.home_city_id - 1]);
            // get all possible paths (up to 2 villages) within distance capacity from home city
            auto allPaths = planAllPaths(data.cities[heli.home_city_id - 1], nearby_villages, heli.distance_capacity);
            // select trips within distance capacity
            auto selectedTrips = selectTripsWithinLimit(allPaths, heli.distance_capacity);
            
            // assign packages to each trip
            for (auto& trip : selectedTrips) {
                trip = assignPackages(trip, heli.weight_capacity, data.packages);
            }
            
            hplan.trips = selectedTrips;
            hplans.push_back(hplan);


        }
        HCState state(hplans);
        return state;
    }

    void add_to_lm(HCState state) {
        b_states.push_back(state);
    }

    HCState estimated_global_extrema() {
        HCState global_state;
        float obj_global = std::numeric_limits<float>::lowest();
        for (auto state : b_states) {
            if (cmp_states(state, global_state, eval_state, red)) {
                global_state = state;
            }
        }
        return global_state;
    }
};



/**
 * @brief to any node, you have to define the state space and the succesor function on your own
 */
class HCState {
public:
    vector<HelicopterPlan> h;

    HCState() {
        ;
    }
    HCState(vector<HelicopterPlan> h) {
        this->h = h;
    }

    vector<HCState> get_successors() { // TODO
        vector<HCState> successors;
        for(HelicopterPlan heli : this->h){
            for(Trip t : heli.trips){
                int length = t.drops.size();
                for(int i = 0; i < length; i++){
                    for(int j = i+1; j < length; j++){
                        swap(t.drops[i], t.drops[j]);
                        successors.push_back(*this);
                    }
                }
            }
        }
        return successors; 
    }

};

/**
 * @brief compares two states and returns true if state1 >= state2
 */
bool cmp_states(HCState state1, HCState state2, float (*func) (HCState), bool red=false) {
    if (!red) {
        if (eval_state(state1) >= eval_state(state2)) {
            return true;
        }
        return false;
    }
    else {
        if (eval_state(state1) <= eval_state(state2)) {
            return true;
        }
        return false;
    }
    return false;
}

/**
 * @brief gets the next succesor (best)
 */
HCState get_best_successor(HCState state, bool red=false) {
    auto successors = state.get_successors();
    float best_obj_value = std::numeric_limits<float>::min();
    HCState best_state = state;
    for (auto successor: successors) {
        float curr_obj_fn = eval_state(successor);
        if (cmp_states(best_state, state, eval_state, red)) {
            best_obj_value = curr_obj_fn;
            best_state = successor;
        }
        
    }
    return best_state;
}

/**
 * @brief produces the objective function for the given state s
 * 
 * @param hcs hill climbing state
 */
float eval_state(HCState hcs) {
    return plan_cost(hcs.h);
}
/**
 * @brief general code for hill climbing with random restarts
 * 
 * @param cstate current state
 * @param bs_state best succsor state
 * @param space search space
 * @param timer termination condition
 * 
 * @remark this restarts under two condition 1. restart eps 2. local extrema
 */
void hcrr(Timer timer, HCState cstate, HCSpace space, bool red=false) {
    if (timer.check_term()) {
        float obj_fn = std::numeric_limits<float>::lowest();  // todo: have to change
        
        HCState bs_state = get_best_successor(cstate);

        // moving to a new search (found a local maxima)
        if (!cmp_states(bs_state, cstate, eval_state, red)) {
            space.add_to_lm(cstate);
            hcrr(timer, space.sample(), space);
        }

        else if (timer.restart()) {
            hcrr(timer, space.sample(), space);
        }

        else{
            hcrr(timer, bs_state,space, red);
        }
    }
    return;
}


/**
 * @brief The main function to implement your search/optimization algorithm.
 * * This is a placeholder implementation. It creates a simple, likely invalid,
 * plan to demonstrate how to build the Solution object. 
 * * TODO: REPLACE THIS ENTIRE FUNCTION WITH YOUR ALGORITHM.
 */
Solution solve(const ProblemData& problem) {
    cout << "Starting solver..." << endl;

    Solution solution;

    // --- START OF PLACEHOLDER LOGIC ---
    // This is a naive example: send each helicopter on one trip to the first village.
    // This will definitely violate constraints but shows the structure.
    
    // timer definition
    float eps_restart = 60;
    Timer timer(problem.time_limit_minutes * 60, eps_restart);

    // maximization problem 
    bool red = false;

    // Seach Space
    HCSpace hcspace(problem, red);

    // Initial Node
    HCState cstate = hcspace.sample();

    // start hill climbing with random restarts
    hcrr(timer, cstate, hcspace, red);

    HCState best_sol = hcspace.estimated_global_extrema();

    for (const auto& helicopter : best_sol.h) {
        HelicopterPlan plan;
        solution.push_back(plan);
    }
    
    // --- END OF PLACEHOLDER LOGIC ---

    cout << "Solver finished." << endl;
    return solution;
}