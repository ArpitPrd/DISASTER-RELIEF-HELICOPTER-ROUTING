#include "solver.h"
#include <iostream>
#include <chrono>

using namespace std;
// You can add any helper functions or classes you need here.

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
float all_trip_cost(Helicopter h) {
    float sum_across_trips = 0;
    for (Trip trip : h.trips) {
        sum_across_trips += trip_cost(trip, h);
    }
    return sum_across_trips;
}

/**
 * @brief general purpose termination criteria
 */
class SW {
public:
    int end_time;
    int ts;
    SW(int end_time) {
        this->end_time = end_time;
        ts = 0;
    }
    bool check_term() {
        return false; // TODO
    }
    float get_time() {
        return 0.0;
    }
    bool restart() {
        return true;
    }
};

/**
 * @brief Defines the search problem space
 * 
 * @param b_states best states
 * @def add_to_lm: add this ""state"" to local maxima
 */
class HCSpace {
public:
    vector<HCState> b_states;
    vector<HCState> p_states;
    HCState sample() {
        ;
    }
    void add_to_lm(HCState state) {
        b_states.push_back(state);
    }
    void add_to_p(HCState state) {
        p_states.push_back(state);
    }
};



/**
 * @brief to any node, you have to define the state space and the succesor function on your own
 */
class HCState {
public:
    vector<Helicopter> h;

    HCState() {
        ;
    }
    HCState(vector<Helicopter> h) {
        this->h = h;
    }

    vector<HCState> get_successors() {
        vector<HCState> successors;
        for(Helicopter heli : this->h){
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

    bool operator==(const HCState state) {
        ;
    }

};

/**
 * @brief compares two states and returns true if state1 >= state2
 */
bool cmp_states(HCState state1, HCState state2, float (*func) (HCState), bool red=true) {
    ;
}

/**
 * @brief gets the next succesor (best)
 */
HCState get_best_successor(HCState state, bool red=true) {
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
 */
float eval_state(HCState h) {
    return 0.0;
}
/**
 * @brief general code for hill climbing with random restarts
 * 
 * @param cstate current state
 * @param bs_state best succsor state
 * @param space search space
 * @param wc termination condition
 */
void hcrr(SW wc, HCState cstate, HCSpace space, bool red=true) {
    if (wc.check_term()) {
        float obj_fn = std::numeric_limits<float>::min();  // todo: have to change
        
        HCState bs_state = get_best_successor(cstate);

        // moving to a new search (found a local maxima)
        if (!cmp_states(bs_state, cstate, eval_state, red)) {
            space.add_to_lm(cstate);
            hcrr(wc, space.sample(), space);
        }

        else if (wc.restart()) {
            space.add_to_p(cstate);
            hcrr(wc, space.sample(), space);
        }

        else{
            hcrr(wc, bs_state,space, red);
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
    
    for (const auto& helicopter : problem.helicopters) {
        HelicopterPlan plan;
        plan.helicopter_id = helicopter.id;

        if (!problem.villages.empty()) {
            Trip trip;
            // Pickup 1 of each package type
            trip.dry_food_pickup = 1;
            trip.perishable_food_pickup = 1;
            trip.other_supplies_pickup = 1;

            // Drop them at the first village
            Drop drop;
            drop.v = problem.villages[0];
            drop.village_id = problem.villages[0].id;
            drop.dry_food = 1;
            drop.perishable_food = 1;
            drop.other_supplies = 1;

            trip.drops.push_back(drop);
            plan.trips.push_back(trip);
        }
        solution.push_back(plan);
    }
    
    // --- END OF PLACEHOLDER LOGIC ---

    cout << "Solver finished." << endl;
    return solution;
}