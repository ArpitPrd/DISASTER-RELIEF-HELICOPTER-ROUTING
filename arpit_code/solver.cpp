#include "solver.h"
#include <iostream>
#include <chrono>
#include <map>
#include <unordered_map>
#include <algorithm>

using namespace std;
using namespace std::chrono;
using hrc = time_point<high_resolution_clock>;
#define now() high_resolution_clock::now()
// You can add any helper functions or classes you need here.

map<int, Helicopter*> hmap;
map<int, Village*> vmap;
map<int, map<int, map<Trip*, bool>>> vill_to_trip_map;
int cnt;



/**
 * @brief calculates the distance in a trip
 */
double trip_distance_travelled(Trip t, Point home) {
    double dist = 0;
    Point prev_coords = home;
    for (Drop d: t.drops) {
        dist += distance(prev_coords, d.v.coords);
        prev_coords = d.v.coords;
    }
    dist += distance(prev_coords, home);
    return dist;
}

float all_trip_distance(vector<Trip> trips, Point home){
    float all_trip_dist = 0;
    for(Trip t : trips){
        float dist = 0;
        Point prev_coords = home;
        for (Drop d : t.drops)
        {
            dist += distance(prev_coords, d.v.coords);
            prev_coords = d.v.coords;
        }
        dist += distance(prev_coords, home);
    }
}

/**
 * @brief use thus for per trip cost
 */
double trip_cost(Trip t, int h_id) {
    if (hmap.find(h_id) == hmap.end()) {
        cerr << "Helicopter ID " << h_id << " not found in hmap!::trip_cost" << endl;
        return 0; // or handle the error as appropriate
    }
    Helicopter h = *hmap[h_id];
    double dist = trip_distance_travelled(t, h.home_city_coords);
    double cost = h.fixed_cost + h.alpha*dist;
    return cost;
}


/**
 * @brief use this for calculation of scores for one trip
 */
float all_trip_cost(HelicopterPlan hplan) {
    if (hmap.find(hplan.helicopter_id) == hmap.end()) {
        cerr << "Helicopter ID " << hplan.helicopter_id << " not found in hmap!::all_trip_cost" << endl;
        return 0; // or handle the error as appropriate
    }
    float sum_across_trips = 0;
    for (Trip trip : hplan.trips) {
        sum_across_trips += trip_cost(trip, hplan.helicopter_id);
    }
    return sum_across_trips;
}

float all_trip_value(HelicopterPlan hplan){
    if (hmap.find(hplan.helicopter_id) == hmap.end()) {
        cerr << "Helicopter ID " << hplan.helicopter_id << " not found in hmap!::all_trip_value" << endl;
        return 0; // or handle the error as appropriate
    }
    float value_across_trips = 0;
    for(Trip trip : hplan.trips){
        value_across_trips += trip.total_value;
    }
    return value_across_trips;
}

float plan_cost(vector<HelicopterPlan> hps) {
    float cost = 0;
    for (auto hplan: hps) {
        cost += all_trip_cost(hplan);
    }
    return cost;
}

float plan_value(vector<HelicopterPlan> hplans){
    float value = 0;
    for(HelicopterPlan hplan : hplans){
        value += all_trip_value(hplan);
    }
    return value;
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
    double end_time;
    hrc last_restart_ts;
    double eps_restart;
    Timer(int end_time, int eps_restart) {
        this->end_time = end_time;
        this->eps_restart = eps_restart;
    }
    void start() {
        start_ts = now();
        last_restart_ts = start_ts;
    }
    bool check_term() {
        hrc ts = now();
        auto duration = duration_cast<seconds>(ts-start_ts);
        // cout << duration.count() << " " << end_time << endl;
        if (duration.count()+40 >= end_time) {
            return true;
        }
        return false; 
    }
    double get_time() {
        hrc ts = now();
        auto duration = duration_cast<seconds>(ts-start_ts);
        return duration.count();
    }
    bool restart() {
        hrc ts = now();
        last_restart_ts = ts;
        auto duration = duration_cast<seconds>(ts-last_restart_ts);
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
 * @param p village from which we want to find other villages within distance dcap
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
        double totalDist = distance(p, v.coords) * 2;
        if (totalDist <= d) {
            // Initialize drop
            Drop drop;
            drop.village_id = v.id;
            drop.v = v;
            drop.dry_food = 0;
            drop.perishable_food = 0;
            drop.other_supplies = 0;
            
            // Initialize trip
            Trip trip;
            trip.dry_food_pickup = 0;
            trip.perishable_food_pickup = 0;
            trip.other_supplies_pickup = 0;
            trip.drops.push_back(drop);
            trip.trip_distance = totalDist;

            // Add to all paths
            allPaths.push_back(trip);
        }
    }

    /*
        Currently removed two-village, since we can construct all from 1 path village
    */

    // // Two-village paths
    // for (size_t i = 0; i < villages.size(); i++) {
    //     for (size_t j = i + 1; j < villages.size(); j++) {
    //         double totalDist = distance(p, villages[i].coords)
    //                          + distance(villages[i].coords, villages[j].coords)
    //                          + distance(villages[j].coords, p);
    //         if (totalDist <= d) {
    //             Drop drop1;
    //             drop1.village_id = villages[i].id;
    //             drop1.v = villages[i];
    //             drop1.dry_food = 0;
    //             drop1.perishable_food = 0;
    //             drop1.other_supplies = 0;   
    //             Drop drop2;
    //             drop2.village_id = villages[j].id;
    //             drop2.v = villages[j];
    //             drop2.dry_food = 0;
    //             drop2.perishable_food = 0;
    //             drop2.other_supplies = 0;
    //             Trip trip;
    //             trip.dry_food_pickup = 0;
    //             trip.perishable_food_pickup = 0;
    //             trip.other_supplies_pickup = 0;
    //             trip.drops.push_back(drop1);
    //             trip.drops.push_back(drop2);
    //             trip.trip_distance = totalDist;
    //             allPaths.push_back(trip);
    //         }
    //     }
    // }

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
Trip assignPackages(Trip& t, int wcap, vector<PackageInfo>& pack) {
    int numVillages = (int)t.drops.size();
    int totalAssigned = 0;
    t.total_value = 0;

    // priority order: perishable (1), dry (0), other (2)
    vector<int> order = {1, 0, 2};

    for (int pkgType : order) {
        int remainingCap = wcap - totalAssigned;
        if (remainingCap <= 0) break;

        // find maximum integer X satisfying constraint
        if (pack.size() <= (size_t) pkgType) {
            cerr << "Package type " << pkgType << " not found in pack!::assignPackages" << endl;
            exit(1);
        }
        int X = remainingCap / (numVillages * pack[pkgType].weight);
        if (X <= 0) continue;

        int assignedThisType = 0;

        // Assign uniformly across all villages
        for (auto& d : t.drops){
            int N = 0;
            if (pkgType == 0){
                if(X > d.v.rem_population_food){
                    d.dry_food += d.v.rem_population_food;
                }
                else{
                    d.dry_food += X;
                }
                N = d.dry_food;
                d.v.rem_population_food -= N;
            }
            if (pkgType == 1){
                if (X > d.v.rem_population_food)
                {
                    d.perishable_food += d.v.rem_population_food;
                }
                else
                {
                    d.perishable_food += X;
                } 
                N = d.perishable_food;
                d.v.rem_population_food -= N;
            }
            if (pkgType == 2){
                if (X > d.v.rem_population_other)
                {
                    d.other_supplies += d.v.rem_population_other;
                }
                else
                {
                    d.other_supplies += X;
                }
                N = d.other_supplies;
                d.v.rem_population_other -= N;
            }
            assignedThisType += X * pack[pkgType].weight;
            t.total_value += N * pack[pkgType].value;
        }

        // Adjust if overshoot
        bool zero_weight = false;
        for (int t =0; t<3; t++){
            if (pack[t].weight <= 0) {
                cerr << "Package weight for type " << t << " is non-positive!::assignPackages" << endl;
                zero_weight = true;
            }
        }
        
        while (assignedThisType > remainingCap && !zero_weight) {
            for (auto& d : t.drops) {
                if (assignedThisType <= remainingCap) break;
                if (pkgType == 0 && d.dry_food > 0) {
                    d.dry_food -= 1;
                    assignedThisType -= pack[pkgType].weight;
                    t.total_value -= pack[pkgType].value;
                }
                if (pkgType == 1 && d.perishable_food > 0) {
                    d.perishable_food -= 1;
                    assignedThisType -= pack[pkgType].weight;
                    t.total_value -= pack[pkgType].value;
                }
                if (pkgType == 2 && d.other_supplies > 0) {
                    d.other_supplies -= 1;
                    assignedThisType -= pack[pkgType].weight;
                    t.total_value -= pack[pkgType].value;
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

bool is_village_visited_by_plan(const HelicopterPlan& hp, int village_id) {
    for (const auto& trip : hp.trips) {
        for (const auto& drop : trip.drops) {
            if (drop.village_id == village_id) {
                return true;
            }
        }
    }
    return false;
}

bool can_add_village_to_trip(const Trip& current_trip, const Village& new_village, const Helicopter& h) {
    // Calculate the distance to the new village and back home.
    Point last_village_coords = current_trip.drops.empty() ?
                                h.home_city_coords :
                                current_trip.drops.back().v.coords;
    
    double new_trip_dist = current_trip.trip_distance - distance(last_village_coords, h.home_city_coords) +
                           distance(last_village_coords, new_village.coords) +
                           distance(new_village.coords, h.home_city_coords);

    if (new_trip_dist > h.distance_capacity) {
        return false;
    }

    return true; // The total Dmax check would be done in the main successor function
}

Drop prepare_drop(Village v) {
    Drop d;
    d.dry_food = 0;
    d.other_supplies = 0;
    d.perishable_food = 0;
    d.village_id = v.id;
    d.v = v;
    return d; // TODO
}

void add_village_to_trip(Trip& trip, const Village& new_village) {
    Drop new_drop = prepare_drop(new_village);
    // The trip's distance will be recalculated later by `assign_packages_to_trip`, so we just add the drop.
    trip.drops.push_back(new_drop);
}

void remove_village_from_trip(Trip& trip, size_t index) {
    if (index < trip.drops.size()) {
        trip.drops.erase(trip.drops.begin() + index);
    }
}

void swap_villages_in_trip(Trip& trip, size_t index1, size_t index2) {
    if (index1 < trip.drops.size() && index2 < trip.drops.size()) {
        std::swap(trip.drops[index1], trip.drops[index2]);
    }
}

bool is_valid_village(int heli_id, int vill_id, Trip& t){
    if(vill_to_trip_map[vill_id][heli_id][&t]){
        return true;
    }
    return false;
}


bool check_dcap(double dcap, Trip t, Point home) {
    double trip_dist = trip_distance_travelled(t, home);
    return dcap >= trip_dist;
}

bool check_dmax(double dmax, vector<Trip> trips, Point home) {
    double total_trips_dist = all_trip_distance(trips, home);
    return total_trips_dist <= dmax;
}



/**
 * @brief to any node, you have to define the state space and the succesor function on your own
 * 
 * @param data problem data originally fed in
 * @note TODO : add a map for each helicopter for each trip maintaining cities travelled in the trip
 */
class HCState {
public:
    vector<HelicopterPlan> hplans;
    ProblemData data;
    
    HCState() {
        ;
    }
    HCState(vector<HelicopterPlan> h, ProblemData data) {
        this->data = data;
        this->hplans = h;
        for(HelicopterPlan heli : h){
            for(Trip t : heli.trips){
                for(Drop d : t.drops){
                    t.vis_villages[d.village_id] = true;
                }
            }
        }
    }

    

    vector<HCState> get_successor() {

        vector<HCState> succ;
        for (size_t h_idx = 0; h_idx < this->hplans.size(); h_idx++) {
            vector<Trip> all_trips = this->hplans[h_idx].trips;
            for (size_t t_idx = 0; t_idx <  all_trips.size(); t_idx++) {
                // we either add a village or remove a village
                Trip expand_this_trip = all_trips[t_idx];
                // adding a village

                vector<Trip> meta_trips = try_adding_village(expand_this_trip, this->hplans[h_idx].helicopter_id, data.d_max);
                for (Trip sug_trip: meta_trips) {
                    HCState hcs_temp(this->hplans, data);
                    assignPackages(sug_trip, hmap[hplans[h_idx].helicopter_id]->weight_capacity, data.packages);
                    hcs_temp.hplans[h_idx].trips[t_idx] = sug_trip;
                    succ.push_back(hcs_temp);
                }
            }

            for(pair<int, Village*> p : vmap){

            }
        }
    }

    vector<HCState> get_successors() { // TODO

        vector<HCState> successors;
        for (size_t h_idx = 0; h_idx < this->hplans.size(); h_idx++) {
            HelicopterPlan heli = this->hplans[h_idx];
            float total = all_trip_cost(heli);
            vector<int> size_check(heli.trips.size(), 0);
            for(size_t t_idx = 0; t_idx < heli.trips.size(); t_idx++){
                Trip t = heli.trips[t_idx];
                Trip trip_temp = t;
                int length = t.drops.size();
                float d = distance(t.drops[length - 1].v.coords, hmap[heli.helicopter_id]->home_city_coords);
                // float cur_trip_cost = trip_cost(t, *hmap[heli.helicopter_id]);
                for(size_t v_idx = 0; v_idx < data.villages.size(); v_idx++){
                    if (trip_temp.vis_villages[data.villages[v_idx].id])
                    {
                        continue;
                    }
                    float d1 = distance(data.villages[v_idx].coords, trip_temp.drops[length - 1].v.coords);
                    float d2 = distance(data.villages[v_idx].coords, hmap[heli.helicopter_id]->home_city_coords);
                    if (((trip_temp.trip_distance - d + d1 + d2) <= hmap[heli.helicopter_id]->distance_capacity) &&
                        (total - d + d1 + d2) <= data.d_max)
                    {
                        Drop new_village;
                        new_village.v = data.villages[v_idx];
                        new_village.village_id = data.villages[v_idx].id;
                        trip_temp.drops.push_back(new_village);
                        trip_temp.vis_villages[data.villages[v_idx].id] = true;
                        trip_temp = assignPackages(t, hmap[heli.helicopter_id]->weight_capacity, data.packages);

                        HCState hcs_temp;
                        hcs_temp.hplans = this->hplans;
                        hcs_temp.hplans[h_idx].trips[t_idx] = trip_temp;
                        hcs_temp.data = data;

                        successors.push_back(hcs_temp);
                    }
                }
                

                // for(int i = 0; i < length; i++){
                //     for(int j = i+1; j < length; j++){
                //         swap(t.drops[i], t.drops[j]);
                //         float new_trip_cost = trip_cost(t, *hmap[heli.helicopter_id]);
                //         if (new_trip_cost <= hmap[heli.helicopter_id]->distance_capacity &&
                //             total - cur_trip_cost + new_trip_cost <= data.d_max){
                //             successors.push_back(*this);
                //         }
                //         swap(t.drops[i], t.drops[j]);
                            
                //     }
                // }

                if(trip_temp.drops.size() == 1){
                    size_check[t_idx] = 1;
                    continue;
                }
                // t.drops.pop_back();
                // t = assignPackages(t, hmap[heli.helicopter_id]->weight_capacity, data.packages);
                // successors.push_back(*this);
                // t = temp;
                
            }

            // creating a new trip
            for(size_t i = 0; i < data.villages.size(); i++){
                Drop d;
                d.v = data.villages[i];
                d.village_id = data.villages[i].id;
                Trip t;
                t.vis_villages[data.villages[i].id] = true;
                t.drops.push_back(d);
                t = assignPackages(t, hmap[heli.helicopter_id]->weight_capacity, data.packages);
                
                HCState hcs_temp;
                hcs_temp.hplans = this->hplans;
                hcs_temp.hplans[h_idx].trips.push_back(t);
                hcs_temp.data = data;

                successors.push_back(hcs_temp);
            }

            // deletion of a trip
            for(size_t i = 0; i < heli.trips.size(); i++){
                if(size_check[i] == 1){
                    HCState hcs_temp;
                    hcs_temp.hplans = hplans;
                    hcs_temp.hplans[h_idx].trips.erase(hcs_temp.hplans[h_idx].trips.begin()+i);
                    hcs_temp.data = data;

                    successors.push_back(hcs_temp);
                }
            }

        }
        return successors; 
    }


};

/**
 * @brief produces the objective function for the given state s
 * 
 * @param hcs hill climbing state
 */
double eval_state(HCState hcs) {
    return plan_value(hcs.hplans) - plan_cost(hcs.hplans);
}

/**
 * @brief compares two states and returns true if state1 >= state2
 */
bool cmp_states(HCState state1, HCState state2, double (*eval_func)(HCState), bool red = false) {
    if (red) {
        return eval_func(state1) <= eval_func(state2);
    }
    return eval_func(state1) >= eval_func(state2);
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
            if (heli.home_city_id <= 0 || (size_t) heli.home_city_id > data.cities.size()) {
                cerr << "City ID " << heli.home_city_id << " not found in data.cities!::sample" << endl;
                exit(1);
            }
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
        HCState state(hplans, data);
        return state;
    }

    void add_to_lm(HCState state) {
        b_states.push_back(state);
    }

    HCState estimated_global_extrema() {
        HCState global_state=b_states[0];
        for (auto state : b_states) {
            if (cmp_states(state, global_state, eval_state, red)) {
                global_state = state;
            }
        }
        return global_state;
    }
};







/**
 * @brief gets the next best successor from a given state
 */
HCState get_best_successor(HCState state, bool red = false) {
    auto successors = state.get_successors();

    // If there are no successors, return the current state
    if (successors.empty()) {
        return state;
    }

    // Initialize best_state with the first successor
    HCState best_state = successors[0];

    // Iterate through the rest of the successors to find the best one
    for (size_t i = 1; i < successors.size(); ++i) {
        if (cmp_states(successors[i], best_state, eval_state, red)) {
            best_state = successors[i];
        }
    }
    
    // An optional improvement is to check if the best successor is actually better than the original state
    if (cmp_states(best_state, state, eval_state, red)) {
        return best_state;
    }
    cout << "No better successor found, return original state" << endl;
    return state; // No better successor found, return original state
}


/**
 * @brief general code for hill climbing with random restarts
 * * @param cstate current state
 * @param space search space
 * @param timer termination condition
 * @param red flag for reduction (minimization)
 */
void hcrr(Timer& timer, HCState cstate, HCSpace& space, bool red = false) {
    // Start the timer
    timer.start();
    
    // Main loop for hill climbing with random restarts
    while (!timer.check_term()) {
        
        // Find the best successor
        HCState bs_state = get_best_successor(cstate, red);

        // Check for a local maximum (no improvement)
        if (!cmp_states(bs_state, cstate, eval_state, red)) {
            space.add_to_lm(cstate); // Add the local maximum to the best states
            cstate = space.sample(); // Random restart
            // cout << "no imp" << endl;
        }
        if (cmp_states(bs_state, cstate, eval_state, red) && cmp_states(cstate, bs_state, eval_state, red)) {
            cout << "same" << endl;
        }
        
        // Check if it's time for a random restart based on interval
        else if (timer.restart()) {
            cstate = space.sample(); // Random restart
        }
        
        // Otherwise, move to the best successor
        else {
            cnt += 1;
            cstate = bs_state; // Climb the hill
            // cout << eval_state(cstate) << endl;
        }
    }

    // After the loop terminates, add the final state to the list of local maxima
    space.add_to_lm(cstate);
}


/**
 * @brief The main function to implement your search/optimization algorithm.
 * * This is a placeholder implementation. It creates a simple, likely invalid,
 * plan to demonstrate how to build the Solution object. 
 * * TODO: REPLACE THIS ENTIRE FUNCTION WITH YOUR ALGORITHM.
 */
Solution solve(const ProblemData& problem) {
    // populating hmap and vmap

    
    for (const auto& h : problem.helicopters) {
        hmap[h.id] = new Helicopter(h);
    }

    for (const auto& v : problem.villages) {
        vmap[v.id] = new Village(v);
    }
    
    
    cout << "Starting solver..." << endl;
    cnt = 0;
    Solution solution;

    // --- START OF PLACEHOLDER LOGIC ---
    // This is a naive example: send each helicopter on one trip to the first village.
    // This will definitely violate constraints but shows the structure.
    
    // timer definition
    double eps_restart = 60;
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
    cout << eval_state(best_sol) << endl;
    solution = best_sol.hplans;
    // --- END OF PLACEHOLDER LOGIC ---

    cout << "Solver finished." << endl;
    return solution;
}