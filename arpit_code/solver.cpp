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

map<int, Village> vmap;
int ascnt = 0;
/**
 * @brief once used information about the packs
 * 
 * @param weight weight
 * @param value value
 */
vector<PackageInfo> packs;


class TripState {
public:
    vector<Drop> drops;
    double wcap;
    double weight;
    TripState(vector<Drop> drops, double wcap, double weight) {
        this->drops = drops;
        this->wcap = wcap;
        this->weight = weight;
    }

    vector<TripState> get_successor() {
        vector<TripState> succs;
        for (size_t d_idx=0; d_idx<drops.size(); d_idx++) {
            for (PackageInfo pack: packs) {
                if (weight + pack.weight <= wcap) {
                    TripState ts(this->drops, this->wcap, this->weight);
                    ts.weight += pack.weight;
                    if (pack.id==0) ts.drops[d_idx].dry_food++;
                    if (pack.id==1) ts.drops[d_idx].perishable_food++;
                    if (pack.id==2) ts.drops[d_idx].other_supplies++;
                    succs.push_back(ts);
                }

                if ((weight-pack.weight>=0) ){
                    TripState ts = *this;
                    ts.weight -= pack.weight;
                    if (pack.id==0 && (ts.drops[d_idx].dry_food>0)){
                        ts.drops[d_idx].dry_food--;
                        succs.push_back(ts);
                    } 
                    if (pack.id==1 && (ts.drops[d_idx].perishable_food>0)){
                        ts.drops[d_idx].perishable_food--;
                        succs.push_back(ts);
                    } 
                    if (pack.id==2 && (ts.drops[d_idx].other_supplies>0)){
                        ts.drops[d_idx].other_supplies--;
                        succs.push_back(ts);
                    }
                    
                }
            }
        }
        
        return succs;
    }

    double eval_state() {
        double value = 0;
        for (Drop d: drops) {
            value += d.dry_food * packs[0].value + d.other_supplies * packs[2].value + d.perishable_food * packs[1].value;
        }
        return value;
    }

    TripState get_best_successor() {
        vector<TripState> succs = get_successor();
        
        TripState best_succ = succs[0];
        
        for (TripState succ: succs) {
            if (succ.eval_state() >= best_succ.eval_state()) {
                best_succ = succ;
            }
        }

        return best_succ;
    }

};

TripState sampler(vector<Drop> drops, double wcap){
    for (Drop &d : drops)
    {
        d.dry_food = 0;
        d.other_supplies = 0;
        d.perishable_food = 0;
    }
    TripState ts(drops, wcap, 0);
    vector<double> pack_density;
    for (const auto pack : packs)
    {
        if (pack.weight > 0)
        {
            pack_density.push_back(pack.density);
        }
        else
        {
            pack_density.push_back(1e9); 
        }
    }
    discrete_distribution<> pack_dist(pack_density.begin(), pack_density.end());

    for(size_t d_idx = 0; d_idx < ts.drops.size(); d_idx++){
        int pack_idx = pack_dist(gen);
        const PackageInfo &pack = packs[pack_idx];

        int population = ts.drops[d_idx].v.population;
        int max_items = 0;
        if(pack_idx == 0 || pack_idx == 1){
            max_items = min(9*population, int((wcap-ts.weight)/pack.weight));
        }
        else{
            max_items = min(population, int((wcap - ts.weight) / pack.weight));
        }
        if(max_items <= 0){
            break;
        }

        uniform_int_distribution<> distribution(0, max_items);
        int assigned_items = distribution(gen);
        ts.weight += assigned_items*pack.weight;
        if(pack_idx == 0){
            ts.drops[d_idx].dry_food += assigned_items;
        }
        if(pack_idx == 1){
            ts.drops[d_idx].perishable_food += assigned_items;
        }
        if(pack_idx == 2){
            ts.drops[d_idx].other_supplies += assigned_items;
        }
    }

    return ts;
}

TripState hcrr_sub(TripState currts, TripState global_state, int rem_restarts, int depth = 0, int max_depth = 1000)
{
    if (depth >= max_depth)
    {
        if(rem_restarts > 0){
            TripState ts = sampler(currts.drops, currts.wcap);
            return hcrr_sub(ts, global_state, rem_restarts - 1, 0, max_depth);
        }
        return global_state;
    }

    TripState best = currts.get_best_successor();

    if (currts.eval_state() >= best.eval_state())
    {
        if(global_state.eval_state() <= currts.eval_state()){
            global_state = currts;
        }
        if(rem_restarts > 0){
            TripState ts = sampler(currts.drops, currts.wcap);
            return hcrr_sub(ts, global_state, rem_restarts-1, 0, max_depth);
        }
        return currts;
    }

    if (global_state.eval_state() <= best.eval_state())
    {
        global_state = currts;
    }

    return hcrr_sub(best, global_state, rem_restarts, depth+1, max_depth);
}

void assignPackages2(Trip &t, double wcap) {
    for (Drop &d: t.drops) {
        d.dry_food = 0;
        d.other_supplies = 0;
        d.perishable_food = 0;
    }
    TripState ts = sampler(t.drops, wcap);
    TripState best = hcrr_sub(ts, ts, 5, 0, 500);
    t.drops = best.drops;

    t.dry_food_pickup = 0;
    t.other_supplies_pickup= 0;
    t.perishable_food_pickup = 0;
    for (Drop &d: t.drops) {
        t.dry_food_pickup += d.dry_food;
        t.perishable_food_pickup += d.perishable_food;
        t.other_supplies_pickup += d.other_supplies;
    }
}



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
    int w = (int) wcap;
  
    vector<PackageInfo> _packs(packs);

    sort(_packs.begin(), _packs.end());

    double total_weight_used = 0;

    vector<int> pops;
    for (Drop &d: drops) {
       
        pops.push_back(d.v.population);
        d.dry_food = 0;
        d.perishable_food = 0;
        d.other_supplies = 0;

    }
    bool improvement_made = true;
    while (improvement_made) {
        improvement_made = false;

        // Iterate through the components, starting with the best one
        for (const auto& pack : _packs) {
            // Find a tuple that can accept this component
            for (size_t d_idx=0; d_idx<drops.size(); d_idx++) {

                bool global_fits = (total_weight_used + pack.weight) <= w;
                if (!global_fits) break; // No need to check local if global fails

                bool local_fits = false;
                if (pack.id == 0 || pack.id == 1) { // Component is a1 or a2
                    if (drops[d_idx].dry_food + drops[d_idx].perishable_food + 1 <= 9 * pops[d_idx]) {
                        local_fits = true;
                    }
                } else { // Component is a3
                    if (drops[d_idx].other_supplies + 1 <= pops[d_idx]) {
                        local_fits = true;
                    }
                }
                if (local_fits) {
                    if (pack.id == 0) drops[d_idx].dry_food+=1;
                    else if (pack.id == 1) drops[d_idx].perishable_food+=1;
                    else drops[d_idx].other_supplies+=1;
                    total_weight_used += pack.weight;
                    improvement_made = true;
                    goto next_iteration; // Restart with the best component
                }
            }
        }
        next_iteration:;
    }

    t.dry_food_pickup = 0;
    t.perishable_food_pickup = 0;
    t.other_supplies_pickup = 0;
    for (Drop &drop: drops) {
        t.dry_food_pickup += drop.dry_food;
        t.perishable_food_pickup += drop.perishable_food;
        t.other_supplies_pickup += drop.other_supplies;
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
 * @note data was deprecated because it is not used in this version
 */
class HCState {
public:
    vector<HelicopterPlan> hplans;

    /*left for empty constructor*/
    HCState() {
        ;
    }

    /**
     * @brief convienience constructor to populate its attrs
     * 
     * @param h provide a vector of all the helicopter plans
     */
    HCState(vector<HelicopterPlan> h) {
        this->hplans = h;
    }

    /**
     * @brief this function calculate the value added across all the possible hps
     * 
     * @return total value added across all the plans
     * @note this must be done on a state level and cannot be done on a trip level, hence here
     */
    double state_value() {
        map<int, Village> vil_map;
        vector<double> food_delivered(vmap.size() + 1, 0.0);
        vector<double> other_delivered(vmap.size() + 1, 0.0);
        vector<double> village_values(vmap.size() + 1, 0.0);
        for (HelicopterPlan hp: hplans) {
            for (Trip t: hp.heli.trips) {
                for (Drop d: t.drops) {
                    int village_id = d.village_id;
                    int vp = d.perishable_food;
                    int vd = d.dry_food;
                    int vo = d.other_supplies;
                    Village village = vmap[d.village_id];

                    double max_food_needed = village.population * 9.0;
                    double food_room_left = max(0.0, max_food_needed - food_delivered[village_id]);
                    double food_in_this_drop = vd + vp;
                    double effective_food_this_drop = min(food_in_this_drop, food_room_left);
                    double effective_vp = min((double)vp, effective_food_this_drop);
                    double value_from_p = effective_vp * packs[1].value;
                    double remaining_effective_food = effective_food_this_drop - effective_vp;
                    double effective_vd = min((double)vd, remaining_effective_food);
                    double value_from_d = effective_vd * packs[0].value;
                    village_values[village_id] += value_from_p + value_from_d;

                    double max_other_needed = village.population * 1.0;
                    double other_room_left = max(0.0, max_other_needed - other_delivered[village_id]);
                    double effective_vo = min((double)vo, other_room_left);
                    village_values[village_id] += effective_vo * packs[2].value;

                    food_delivered[village_id] += food_in_this_drop;
                    other_delivered[village_id] += vo;
                }
            }
        }

        double value = accumulate(village_values.begin(), village_values.end(), 0.0);

        return value;
    }

    /**
     * @brief calculates the cost for my state
     * 
     * @return cost of my state in double
     */
    double state_cost() {
        double cost = 0;
        for (HelicopterPlan hp: hplans) {
            cost += hp.heli.trips_cost(hp.heli.trips);
        }
        return cost;
    }

    /**
     * @brief calculate the objective function of my state by value-cost
     * 
     * @return the objective function of my state
     */
    double eval_state() {
        return state_value() - state_cost();
    }
    

    
    /**
     * @brief finds the successors or the children only, this does not contain me
     * 
     * @return a vector of HCState that are my children
     */
    vector<HCState> get_successor() {

        vector<HCState> succ;
        for (size_t h_idx = 0; h_idx < hplans.size(); h_idx++) {
            vector<Trip> all_trips = hplans[h_idx].heli.trips;
            vector<bool> trip_size_one(all_trips.size());
            
            for (size_t t_idx = 0; t_idx <  all_trips.size(); t_idx++) {
                
                // adding a village
                Trip t = all_trips[t_idx];
                
                vector<Trip> meta_trips = hplans[h_idx].heli.try_adding_village(t_idx, vmap);
                
                for (Trip sug_trip: meta_trips) {
                    HCState hcs_temp(hplans);
                    assignPackages2(sug_trip, hplans[h_idx].heli.weight_capacity);
                    hcs_temp.hplans[h_idx].heli.trips[t_idx] = sug_trip;
                    succ.push_back(hcs_temp);
                }

                Trip rem_vill = hplans[h_idx].heli.try_removing_village(t_idx, vmap);
                if(rem_vill.drops.size() != 0){
                    HCState hcs_temp(hplans);
                    assignPackages2(rem_vill, hplans[h_idx].heli.weight_capacity);
                    hcs_temp.hplans[h_idx].heli.trips[t_idx] = rem_vill;
                    succ.push_back(hcs_temp);
                }
                else{
                    trip_size_one[t_idx] = true;
                }
 
            }

            vector<Trip> one_village_trips = hplans[h_idx].heli.try_new_trip(vmap);

            for (Trip sug_trip: one_village_trips) {
                HCState hcs_temp(hplans);
                assignPackages2(sug_trip, hplans[h_idx].heli.weight_capacity);
                hcs_temp.hplans[h_idx].heli.trips.push_back(sug_trip);
                succ.push_back(hcs_temp);
            }

            for (size_t t_idx = 0; t_idx < all_trips.size(); t_idx++){
                Trip t = all_trips[t_idx];
                if(trip_size_one[t_idx]){
                    HCState hcs_temp(hplans);
                    all_trips.erase(all_trips.begin() + t_idx);
                    hcs_temp.hplans[h_idx].heli.trips = all_trips;
                    succ.push_back(hcs_temp);
                    all_trips.insert(all_trips.begin() + t_idx, t);

                }
            }
        }
        return succ;
    }


    /**
     * @brief compares this state with someother state and returns true if I am better
     * 
     * @param state to compare me with
     * @return true if I am better
     */
    bool cmp_states(HCState state) {
        return this->eval_state() >= state.eval_state();
    }
    
    /**
     * @brief wraps around get_successor to find the best child
     * 
     * @return HCState type of the next succesor you can potentially explore
     * 
     * @warning does not compare myself with the best that is found
     */
    HCState get_best_successor() {
        vector<HCState> succs = get_successor();
        if(succs.size() == 0){
            
        }
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
 * 
 * The following helper functions are defined
 * 
 * villagesWithinDistance finds all the villages within a particular distance about a particyular center
 * prepare_drop a simple template preparer
 * prepare_trip a simple trip preparer
 * sample a function to sample from the search space
 * add_to_lm adds to local maxima vector for convinience
 * estimated_global_extrema finds the best solution so far in terms of the objective function
 * 
 * @remark sampler is put here since you "sample" from a "space"
 * @note this acts as a wrapper to the entire space of solutions
 */
class HCSpace {
public:
    vector<HCState> b_states;
    ProblemData data;
    HCState global;
    
    /**
     * @brief convinience wrapper to populating the parameters of the state
     * 
     * @param _data pass the problem data wrapper
     */
    HCSpace(ProblemData _data) {
        data = _data;
    }

    /**
     * @brief get all villages within distance dcap from center village
     * 
     * @param villages list of all villages
     * @param p village from which we want to find other villages within distance dcap
     * @param dcap distance capacity
     * @return std::vector<Village> list of villages within distance dcap from center village
     * 
     * @note no change done to villages, p
     * @remark this just calculates all the villages within a circle of radious dcap, do not modify any firther
     * @remark trip distance is twice of one way
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

    /**
     * @brief produces a sample state HCState from random sampling technique mentioned below
     * 
     * 0. for each helicopter we find a plan as follows
     * 1. find all the villages within the range of the helicopter (dcap get satisfied)
     * 2. shuffle these villages in a random order
     * 3. collect a subset of villages randomly such that Dmax is satisfied
     * 4. name the trips (vil) for all the vil that are selected
     * 5. wrap all these trips into a Helicopter
     * 6. wrap all the helicopters plans into HCState
     * 7. return this state
     * 
     * @return HCState which is valid
     */
    HCState sample() {
        HCState state;
        for (const auto& heli : data.helicopters) {

            HelicopterPlan hplan;
            hplan.heli = heli;
            hplan.helicopter_id = heli.id;
            hplan.helicopter_id = heli.id;

            vector<Village> v_in_range = villagesWithinDistance(data.villages, heli.distance_capacity/2, heli.home_city_coords);
            
            shuffle(v_in_range.begin(), v_in_range.end(), gen);

            double d_tot = 0;
            for (Village v_pros: v_in_range) {
                d_tot += 2*distance(v_pros.coords, heli.home_city_coords);
                if (d_tot <= heli.d_max) {
                    Trip t = prepare_trip(v_pros);
                    assignPackages2(t, hplan.heli.weight_capacity);
                    hplan.heli.trips.push_back(t);
                }
            }

            state.hplans.push_back(hplan);

        }
        return state;
    }

    /**
     * @brief adds to the local maxima list
     * 
     * @param state local maxima found
     * @note use carefully do not add any random node
     */
    void add_to_lm(HCState state) {
        b_states.push_back(state);
    }

    /**
     * @brief finds the global maxima across all the local maximas found
     * 
     * @return HCState the global maxima
     */
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
        if (duration.count() + 10 >= end_time) {
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
 * 
 * @param cstate current state
 * @param space search space
 * @param timer termination condition
 * @param red flag for reduction (minimization)
 * 
 */
void hcrr(Timer& timer, HCState cstate, HCSpace& space) {
    // Start the timer
    timer.start();
    
    // Main loop for hill climbing with random restarts
    while (!timer.check_term()) {
        
        HCState bs_state = cstate.get_best_successor();
        
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
        }
    }

    // After the loop terminates, add the final state to the list of local maxima
    space.add_to_lm(cstate);
}

/**
 * @brief The main function to implement your search/optimization algorithm.
 * * This is a placeholder implementation. It creates a simple, likely invalid,
 * plan to demonstrate how to build the Solution object. 
 * 
 * @param problem data constructor to the problem
 */
Solution solve(const ProblemData& problem) {
    
    ascnt = 0;
    for (Village v: problem.villages) {
        vmap[v.id] = v;
    }


    packs = problem.packages;
    for (int i=0; i<=2; i++) {
        packs[i].density = packs[i].value / packs[i].weight;
        packs[i].id = i;
    }
    
    
    //cout << "Starting solver..." << endl;
    Solution solution;
    
    // timer definition
    double eps_restart = 60;
    Timer timer(problem.time_limit_minutes * 60, eps_restart);

    // Seach Space
    HCSpace hcspace(problem);

    // Initial Node
    HCState cstate = hcspace.sample();

    // start hill climbing with random restarts
    hcrr(timer, cstate, hcspace);

    HCState best_sol = hcspace.estimated_global_extrema();
    //cout << "score: " << best_sol.eval_state() << endl;
    solution = best_sol.hplans;

    //cout << "Solver finished." << endl;
    return solution;
}