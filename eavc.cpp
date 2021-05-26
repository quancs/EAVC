#include "eavc.h"

static unsigned int eaSelTimes = 24;
void eavc_search() {
    unsigned removeV, addV;
    unsigned e, v1, v2;
    step = 1;

    while (1) {
        if (ueSize == 0) { //update best solution if needed
            update_best_sol();
            //if (c_size==optimal_size) return;
            remove_one_vertex_from_c();
            continue;
        }

        if (step % 16 == 0) { //check cutoff
#ifdef WIN32
            finish = clock();
            long elapsedClocks = finish - gstart;
#else
            times(&finish);
            long elapsedClocks = finish.tms_utime + finish.tms_stime - startTime;
#endif
            if (elapsedClocks >= clocksTotal)
                return;
        }

        removeV = choose_remove_v();
        remove(removeV);

        //EA strategy
        if (ueSize == 1)
            e = uncovEdges[0];
        else {
            e = uncovEdges[static_cast<unsigned int>(pcg32_fast()) % ueSize];
            unsigned int et = eTimestamp[e];
            for (unsigned int i = 0; i < eaSelTimes; i++) {
                unsigned int e2 = uncovEdges[static_cast<unsigned int>(pcg32_fast()) % ueSize];
                unsigned int et2 = eTimestamp[e2];
                if (et2 < et) {
                    e = e2;
                    et = et2;
                }
            }
        }
        //e = uncovEdges[static_cast<unsigned int>(pcg32_fast()) % ueSize];

        v1 = edges[e].v1;
        v2 = edges[e].v2;

        if (dscore[v1] > dscore[v2] || (dscore[v1] == dscore[v2] && vTimestamp[v1] < vTimestamp[v2]) )
            addV = v1;
        else
            addV = v2;

        add(addV);

        unsigned int index = cVertexIndexes[removeV];
        cVertexIndexes[removeV] = 0;

        cVertexes[index] = addV;
        cVertexIndexes[addV] = index;

        vTimestamp[addV] = vTimestamp[removeV] = step;

        //cout << "c Best vertex cover size = " << best_c_size <<  ", SearchSteps = " << best_step <<  ", SearchTime = " << best_comp_time << endl;
        step++;
    }
}



int main(int argc, char *argv[]) {
    unsigned int seed, i;
    if (build_instance(argv[1]) != 1) {
        cout << "can't open instance file" << endl;
        return -1;
    }

    //optimal_size=0;
    i = 2;

    sscanf(argv[i++], "%d", &seed);
    sscanf(argv[i++], "%d", &cutoffTime);
    sscanf(argv[i++], "%d", &eaSelTimes);

    pcg32_fast_init(seed);

    cout << "c EAVC: " << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[4] << endl;

#ifdef WIN32
    gstart = clock();
    clocksTotal = static_cast<int>(cutoffTime) * CLOCKS_PER_SEC;
#else
    times(&gstart);
    startTime = gstart.tms_utime + gstart.tms_stime;
    clocksTotal = cutoffTime * sysconf(_SC_CLK_TCK);
#endif
    preprocess();
    unsigned initSize = 0;
    if (ueSize > 0) {
        init_sol();
        initSize = bestSize;
        eavc_search();
        cVertexes = cAllVertexes;
    } else
        cout << "c no uncovered edge after preprocessing" << endl;
    //check solution
    if (check_solution() == 1) {
        //print_solution();
        //compute times
#ifdef WIN32
        bestTime = bestTime / CLOCKS_PER_SEC;
        finish = clock();
        double totalTime = (static_cast<double>(finish - gstart)) / CLOCKS_PER_SEC;
#else
        bestTime = bestTime / sysconf(_SC_CLK_TCK);
        times(&finish);
        double totalTime = double(finish.tms_utime - gstart.tms_utime + finish.tms_stime - gstart.tms_stime) / sysconf(_SC_CLK_TCK);
#endif
        bestTime = round(bestTime * 100) / 100.0;

        cout << "c steps total:" << step << endl;
        cout << "c time  total:" << totalTime << endl;
        cout << "c Best vertex cover size = " << bestSize << ", init = " << initSize <<  ", steps = " << bestStep <<  ", time = " << bestTime << endl;
    }

    free_memory();
    return 0;
}
