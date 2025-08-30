#include "solver.h"
#include <iostream>
#include <chrono>
#include <map>

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
        ;
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

    vector<HCState> get_successors() {
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