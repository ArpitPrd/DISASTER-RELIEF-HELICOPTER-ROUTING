#include "solver.h"
#include <iostream>
#include <chrono>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <random>


using namespace std;
using namespace std::chrono;
std::random_device rd;
std::mt19937 gen(rd());
using hrc = time_point<high_resolution_clock>;
#define now() high_resolution_clock::now()

/**
 * @brief once used information about the packs
 * 
 * @param weight weight
 * @param value value
 */
vector<PackageInfo> packs;

/**
 * @brief assigns packages by optimising the following
 * 
 * sigma over i v1*ai1 + v2*ai2 + v3*ai3
 * 
 * such that,
 * ai1 + ai2 <= 9*Ni and ai3 <= Ni where i is the ith village in the sequence of all the drops
 * 
 * @param t trip for which we find the drops
 * @param wcap weight capacity of the trip
 * 
 * @return none, modified in place
 * 
 * @note highly expensive function, need to think of alternates
 */
void assignPackages(Trip &t, double wcap) {
    vector<Drop> drops = t.drops;
    int n = t.drops.size();
    int w = (int) wcap;

    vector<PackageInfo> _packs(packs);

    sort(_packs.begin(), _packs.end());

    long long total_weight_used = 0;

    vector<int> pops;
    for (Drop d: drops) {
        pops.push_back(d.v.population);
    }
    bool improvement_made = true;
    while (improvement_made) {
        improvement_made = false;

        // Iterate through the components, starting with the best one
        for (const auto& pack : _packs) {
            // Find a tuple that can accept this component
            for (int d_idx=0; d_idx<drops.size(); d_idx++) {

                Drop drop = drops[d_idx];
                // Check global constraint
                bool global_fits = (total_weight_used + pack.weight) <= w;
                if (!global_fits) continue; // No need to check local if global fails

                // ***MODIFIED: Check the new, specific local constraints***
                bool local_fits = false;
                if (pack.id == 0 || pack.id == 1) { // Component is a1 or a2
                    if (drop.dry_food + drop.perishable_food + 1 <= 9 * pops[d_idx]) {
                        local_fits = true;
                    }
                } else { // Component is a3
                    if (drop.other_supplies + 1 <= pops[d_idx]) {
                        local_fits = true;
                    }
                }

                if (local_fits) {
                    // Add the component to this tuple
                    if (pack.id == 0) drop.dry_food++;
                    else if (pack.id == 1) drop.perishable_food++;
                    else drop.other_supplies++;
                    
                    total_weight_used += pack.weight;
                    improvement_made = true;
                    goto next_iteration; // Restart with the best component
                }
            }
        }
        next_iteration:;
    }

    t.drops = drops;

    return;
 }

/**
 * @brief to any node, you have to define the state space and the succesor function on your own
 * 
 * this is a wrapper to the Solution supporting various functions listed here
 * 
 * state_value - provides the value given of this state
 * state_cost - provides the net cost expended in this state
 * eval_state - provides the numeric good-ness of this function
 * get_successor - provides the list of successor based on simple addition of new villages in a trip
 * get_best_successor - from get_successors provides the best succesor in the children (does not compare with itself also)
 * cmp_state - compare itself with some give state and returns true if I am better than the given
 * 
 * @param hplans all the helicopter plans this associates a state
 * @param data problem data originally fed in
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
    }

    /**
     * @brief this function calculate the value added across all the possible hps
     * 
     * @return total value added across all the plans
     * @note this must be done on a state level and cannot be done on a trip level
     */
    double state_value() {
        map<int, Village> vil_map;

        for (HelicopterPlan hp: hplans) {
            for (Trip t: hp.heli.trips) {
                for (Drop d: t.drops) {
                    
                    // if village found for the first time
                    
                    if (vil_map.find(d.village_id)==vil_map.end()) {
                        vil_map[d.village_id] = d.v;
                    }

                    // if already exists
                    else {
                        vil_map[d.village_id].pack_received.dry_food += d.dry_food;
                        vil_map[d.village_id].pack_received.perishable_food += d.perishable_food;
                        vil_map[d.village_id].pack_received.other_supplies += d.other_supplies;
                    }

                }
            }
        }

        // now eval for each villages
        double value = 0;
        for (pair<int, Village> p: vil_map) {
            value += p.second.value_receieved(packs);
        }

        return value;
    }

    /**
     * uses the inbuilt function to find the cost of each hp
     * 
     * this is suppored per hp
     */
    double state_cost() {
        double cost = 0;
        for (HelicopterPlan hp: hplans) {
            cost += hp.trips_cost();
        }

        return cost;
    }

    double eval_state() {
        return state_value() - state_cost();
    }
    

    

    vector<HCState> get_successor() {

        vector<HCState> succ;
        for (size_t h_idx = 0; h_idx < this->hplans.size(); h_idx++) {
            vector<Trip> all_trips = this->hplans[h_idx].heli.trips;
            for (size_t t_idx = 0; t_idx <  all_trips.size(); t_idx++) {
                
                // adding a village
                Trip t = all_trips[t_idx];
                vector<Trip> meta_trips = hplans[h_idx].heli.try_adding_village(t);
                for (Trip sug_trip: meta_trips) {
                    HCState hcs_temp(this->hplans, data);
                    assignPackages(sug_trip, hplans[h_idx].heli.weight_capacity);
                    hcs_temp.hplans[h_idx].heli.trips[t_idx] = sug_trip;
                    succ.push_back(hcs_temp);
                }
            }
        }
        return succ;
    }


    /**
     * @brief compares this state with someother state and returns true if I am better
     */
    bool cmp_states(HCState state) {
        return this->eval_state() >= state.eval_state();
    }
    
    HCState get_best_successor() {
        vector<HCState> succs = get_successor();

        HCState best_hcs = succs[0];

        for (HCState succ: succs) {
            if (succ.cmp_states(best_hcs)) {
                best_hcs = succ;
            }
        }
        return best_hcs;
    }

};


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

    Drop prepare_drop(Village v) {
        Drop d;
        d.dry_food = 0;
        d.other_supplies = 0;
        d.perishable_food = 0;
        d.v  = v;
    }

    Trip prepare_trip(Village v) {
        Trip t;
        Drop d = prepare_drop(v);
        t.drops = {d};
        t.dry_food_pickup = 0;
        t.other_supplies_pickup = 0;
        t.perishable_food_pickup = 0;
    }

    HCState sample() {
        vector<HelicopterPlan> hplans;
        for (const auto& heli : data.helicopters) {
            HelicopterPlan hplan;
            hplan.helicopter_id = heli.id;

            vector<Village> v_in_range = villagesWithinDistance(data.villages, heli.distance_capacity, heli.home_city_coords);

            shuffle(v_in_range.begin(), v_in_range.end(), gen);

            double d_tot = 0;
            for (Village v_pros: v_in_range) {
                d_tot += distance(v_pros.coords, heli.home_city_coords);
                if (d_tot < heli.d_max) {
                    Trip t = prepare_trip(v_pros);
                    hplan.heli.trips = {t};
                }
            }

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
            if (state.cmp_states(global_state)) {
                global_state = state;
            }
        }
        return global_state;
    }
};


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
        HCState bs_state = cstate.get_best_successor();

        // Check for a local maximum (no improvement)
        if (cstate.cmp_states(bs_state)) {
            space.add_to_lm(cstate); // Add the local maximum to the best states
            cstate = space.sample(); // Random restart
        }
        
        // Check if it's time for a random restart based on interval
        else if (timer.restart()) {
            cstate = space.sample(); // Random restart
        }
        
        // Otherwise, move to the best successor
        else {
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

    packs = problem.packages;
    for (int i=0; i<=2; i++) {
        packs[i].density = packs[i].value / packs[i].weight;
        packs[i].id = i;
    }
    
    
    cout << "Starting solver..." << endl;
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
    cout << best_sol.eval_state() << endl;
    solution = best_sol.hplans;
    // --- END OF PLACEHOLDER LOGIC ---

    cout << "Solver finished." << endl;
    return solution;
}