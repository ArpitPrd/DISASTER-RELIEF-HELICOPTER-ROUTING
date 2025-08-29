## Problem Statement
- Complex Problem and formulate as search 
- Heuristic: Optimal 
- Local: Good for large problems

- Urgent need to carry out relief ops, need neighbouring cities
- Cities to Villages via H for N packages

- Goal: Find the allocation of relief goods to villages under constraints of H, maximising aid delieverd + minimizing "logistical strain cost", which is caused by routing inefficiently

- this forms the glue we want to stretch but are not able to

- Information on villages:
    - (xv, yv)
    - distance: euclidean (city-city or city-village) (km)
    - nv: number of people expected to be stranded

- information on cities:
    - (xc, yc)
    - distance: euclidean (city-city or city-village) (km, dk may be useful later)

- information on packets:
    - dry food
    - wet (perishable food)
    - other supplies
    - w(t): fixed weight of the packet, e.g. weight(wet food) = 2
    - v(t): fixed value of the packet, e.g. value(dry food)=10


H Helicopters
- Information on Helicopters
    - home(h): Home City of the helicopter
    - wcap(h): weight capacity per trip
    - dcap(h): distance capacity per trip

- Quantitative Goal: 
    - send 9 meals per stranded person 
        - dry/ wet works
        - wet > dry preference
    - about 1 unit of other supplies per stranded person 
    - for all Helicopter h:
        - h must start and end at home city (this is the definition of a trip)
        - w <= wcap(h)
        - d <= dcap(h)
        - sigma di <= Dmax
            - di = distance travelled in the ith trip from home to home
    - max (total value achieved - total trip cost)
        - total value achieved = sigma v(all food items delivered)
        - total trip cost = sigma across all trips for each flight and across flights (f + alpha * distance per trip)
    Find:
        - for each Helicopter h:
            print(len h.trip)
            for each trip t in h.trip:
                for each package p in h.trip.package:
                    print p.name, p.freq
                for each village v in InOrder(h.trip.villages):
                    print v.name
                    for each package p in v.package:
                        print p.name, p.freq

- Description of Input
    - total processing time (can vary this to change the alg)
    - Dmax
    - w(d) v(d) w(p) v(p) w(o) v(o)
    - C xc1 yc1 xc2 yc2 xc3 yc3 ... xcC ycC
    - V xv1 yv1 nv1 xv2 yv2 nv2 ... xvV yvV nvV
    - H cityid1 wcap1 dcap1 F1 alpha1 cityid2 wcap2 dcap2 F2 alpha2 ... cityidH wcapH dcapH FH alphaH

    - check the submission details

- Desciprtion of the output
```
    for each Helicopter h:
        print h.id len(h.trip) end="\n"
        for each trip t in h.trip:
            for each package p in h.trip.package:
                print p.freq end=" "
            print len(h.trip.villages) end=" "
            for each village v in InOrder(h.trip.villages):
                print v.id end=""
                for each package p in v.package:
                    print p.freq

        print -1
```
- Evalutation of each cost
```
    cost = 0
    // cost subtraction
    for each Helicopter h:
        sum_across_trips = 0
        for each trip t in h.trip:
            dist_total = 0
            for each village v in InOrder(h.trip.villages):
                dist_total += dist(v.c, prev_c)
                for each packkage p in v.packages:
                    cost += p.freq * p.v // write a separate function
            sum_across_trips + h.F + h.alpha * dist_total
        cost -= sum_across_trips
    
    fucntion eval_cost
    // contains the check if more than 9 are considered, the consider the best ones

    TODO
    can check just the number of excesses, others just check p.v for perisible if not exceed then (this will bring us an intial working code)

    // cost addition
    for each Helicopter h:
        for each trip t in h.trip:
            for each village v in InOrder(h.trip.villages):
                cost += eval_cost(v.packages)
```
- Verification
```
    // distance 
    for each Helicopter h:
        h_trip_dist = 0
        for each trip t in h.trip:
            dist_total = 0
            for each village v in InOrder(h.trip.villages):
                dist_total += dist(v.c, prev_c)
                h_trip_dist += dist(v.c, prev_c)
            assert dist_total <= h.dcap
        assert h_trip_dist <= DMax
        
    // weight
    for each Helicopter h:
        for each trip t in h.trip:
            weight_total = 0
            for each village v in InOrder(h.trip.villages):
                for each package p in v.packages:
                    weight_total += p.weight * p.freq
            assert weight_total <= h.wcap
```

## Algorithms to test

- Heuristic Search: Quality may be good, but very slow
- Branch and Bound: 
    - state space
    - Transition function
    - eval a heuritci fucntion 
    - Depth First Branch and Bound
    - Not scalable
- Local Search (reco)
    - quality of search as a function of time

## Submission

- VMs
    - amd64
    - g++ 9.4.0
    - RAM 8GB

- check the input is valid or not
    1
    100
    0.01 1 0.1 2 0.005 0.1
    2 0 0 10 10
    2 0 5 1000 0 10 1000
    2 1 100 25 10 1 2 100 50 10 1

- submission:
    - zipping:
        - Name: 2022EE11837_2022EE111XX.zip
        - unzipping should produce:
            - (in working dir)
            - .cpp/ .h files
            - Makefile
            - writeup.txt
    - make produces the executable "main"
    - ./main input.txt output.txt

    -writeup.txt
        - Line 1: Even though I/we have taken help from the following
        students and LLMs in terms of discussing ideas and coding practices, all my/our code is written by
        me/us
        - Students:
        - LLMs
        - Something about the code

- 12hr log instructions

- critteria of judgment: objective function, how well it performs

## Search Probelm Formualtion

- Modelling
    - A state here would include:
        - Helicopters: their trips and order of village travels

- Algortithm

## General Pointer

- c++: cannot define recursilvely in same class
- std::next(n, it), gives an iterator n steps forwards'
- bool operator==(const template& var) {}
- can define this inside a struct also
- pass function as an arg as int (*func)(int, int)

## Coding Enhancements

- >= operation not in class, because it requires the obj function, not truly a class property

## Plans

- 28 Aug
    - Input Output formatting V
    - Hill climbing code A
    - Modelling Tomorrow