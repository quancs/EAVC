#ifndef EAVC_H
#define EAVC_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string.h>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <functional>
#include <set>
#include <cassert>
#include "pcg32.h"
using namespace std;

#ifdef WIN32
#include <time.h>
static clock_t start, finish, gstart;
#else
#include <unistd.h>
#include <sys/times.h>
static tms start, finish, gstart;
static int startTime;
#endif
static long clocksTotal;

struct Edge {
    unsigned int v1;
    unsigned int v2;
};

//possible status of vertice
enum VerticeStatus : unsigned char {
    NIN = 0/*not in cover*/, NIN_ABS = 1/*not in cover absolutely*/, UNKNOWN = 2, IN = 3/*in cover*/, IN_ABS = 4 /*in absolutely*/
};

/*parameters of algorithm*/
//static unsigned int	maxSteps; //step limit
static unsigned int	cutoffTime; //time limit
static unsigned int	step;
//static unsigned int	optimalSize; //terminate the algorithm before step limit if it finds a vertex cover of optimal_size

/*structures about graph*/
static unsigned int vNum;//|V|: 0...v-1
static unsigned int eNum;//|E|: 0...e-1
static Edge *edges;
static unsigned int	**vEdges;//edges related to v, v_edges[i][k] means vertex v_i's k_th edge
static unsigned int	**vAdj;//v_adj[v_i][k] = v_j(actually, that is v_i's k_th neighbor)
static unsigned int	*vDegree;//amount of edges (neighbors) related to v

/*structures of preprocessing */
static unsigned *dVDegree;//dynamic vertice degree
static unsigned *uvSet;//untested vertice set
static unsigned uvSize;//untested vertice size
static VerticeStatus *checkAdjDomV;//if >0, it means we has the necessity to check if adj dominates v, checkAdjDomV=vInC
static unsigned total_IN_ABS;//total IN_ABS vertices
static unsigned total_NIN_ABS;//total NIN_ABS vertices

/* structures about solution */
//current candidate solution
static VerticeStatus *vInC;//a status indicates whether a vertex is in C
static unsigned int cAVSize;//cardinality of IN_ABS and IN vertices in cAllVertexes
static unsigned int *cAllVertexes;//a backup of cVertexes, layout:IN_ABS vertices | IN vertices | unused area
static unsigned int cvSize;//cardinality of cVertexes
static unsigned int	*cVertexes;//an array consists of only vertices in C, layout:IN vertices | unused area
static unsigned int	*cVertexIndexes;//if v belongs to C, cVertexes[cVertexIndexes[v]]=v
static int	*dscore;//dscore of v
static unsigned int *vTimestamp;
//uncovered edges of current solution
static unsigned int	*uncovEdges;//store the uncov edge number
static unsigned int ueSize;
static unsigned int	*uncovEdgeIndexes;//which position is an edge in the uncov_stack
static unsigned int *eTimestamp;

//best solution found
static unsigned int	bestSize;
static VerticeStatus *vInBest;//a flag indicates whether a vertex is in best solution
static double bestTime;
static unsigned int bestStep;

/* functions declaration */
int build_instance(char *filename);
/* preprocess functions: after preprocessing, dVDegree contains all the vertice degrees */
void prep_remove_from_neighbours(unsigned v);
bool prep_handle_degree_0(unsigned v);
bool prep_handle_degree_1(unsigned v);
bool prep_handle_degree_2(unsigned v);
bool prep_handle_dom(unsigned v, bool &found);
void preprocess();
/* initiate functions*/
void init_sol();
void init_vertex_greedy_vc(int *dscore1, VerticeStatus *vInC1, unsigned int &cvSize1);
void init_edge_greedy_vc(int *dscore2, VerticeStatus *vInC2, unsigned int &cvSize2);//EdgeVC
void init_match_vc(int *dscore3, VerticeStatus *vInC3, unsigned int &cvSize3);//MatchVC
/* searching related functions*/
void eavc_search();
void add(unsigned int v);
void remove(unsigned int v);
int check_solution();
bool check_solution(VerticeStatus *vInCTest, unsigned int cvSizeTest, int *scoreTest);
bool in(VerticeStatus *vs, unsigned v);//test if v in vs
bool not_in(VerticeStatus *vs, unsigned v);//test if v not in vs


inline void update_best_sol() {
    memcpy(vInBest, vInC, sizeof(VerticeStatus)*vNum);

    bestSize = total_IN_ABS + cvSize;
#ifdef WIN32
    finish = clock();
    bestTime = (static_cast<double>(finish - gstart));
#else
    times(&finish);
    bestTime = double(finish.tms_utime + finish.tms_stime - startTime);
#endif
    bestStep = step;
    //cout << "c Best vertex cover size = " << best_c_size <<  ", SearchSteps = " << best_step <<  ", SearchTime = " << best_comp_time << endl;
}

int build_instance(char *filename) {
    char line[1024];
    char tempstr1[10];
    char tempstr2[10];
    unsigned int  v, e;

    char	tmp;
    unsigned int  v1, v2;

    ifstream infile(filename);
    if (!infile.is_open()) return 0;

    /*** build problem data structures of the instance ***/
    infile.getline(line, 1024);
    while (line[0] == '%') infile.getline(line, 1024);
    if (line[0] == 'p')
        sscanf(line, "%s %s %d %d", tempstr1, tempstr2, &vNum, &eNum);
    else
        sscanf(line, "%d %d %d", &v, &vNum, &eNum);

    eTimestamp = new unsigned int[eNum];
    edges = new Edge [eNum];							//be initialized here
    uncovEdges = new unsigned int[eNum];                      //only need to initialized uncov_stack_fill_pointer, has been done in init_sol()
    uncovEdgeIndexes = new unsigned int [eNum];             //the same as above
    dscore = new int [vNum];                       //be initialized in init_sol()
    vTimestamp = new unsigned int [vNum];             //be initialized in init_sol()
    vEdges = new unsigned int *[vNum];                     //be initialized here
    vAdj = new unsigned int *[vNum];                       //the same as above
    vDegree = new unsigned int [vNum];                     //the same as above
    memset(vDegree, 0, sizeof(unsigned int) * vNum);
    vInC = new VerticeStatus [vNum];                      //be initialized in init_sol()
    cVertexes = new unsigned int [vNum];                  //be initialized in reset_remove_cand() in init_sol()
    cVertexIndexes = new unsigned int [vNum];         //the same as above
    vInBest = new VerticeStatus [vNum];                 //be initialized in update_best_sol() in init_sol()

    for (e = 0; e < eNum; e++) {
        infile.getline(line, 1024);
        stringstream ss;
        ss << line;
        if (line[0] == 'e')
            ss >> tmp >> v1 >> v2;
        else
            ss >> v1 >> v2;

        v1 = v1 - 1;
        v2 = v2 - 1;

        vDegree[v1]++;
        vDegree[v2]++;

        edges[e].v1 = v1;
        edges[e].v2 = v2;
        eTimestamp[e] = 0;
    }
    infile.close();

    /* build v_adj and v_edges arrays */
    for (v = 0; v < vNum; v++) {
        vAdj[v] = new unsigned int[vDegree[v]];
        vEdges[v] = new unsigned int[vDegree[v]];
    }

    //for(v=1; v<=v_num; v++) v_degree_tmp[v]=0;

    //	int v_degree_tmp[MAXV];
    int *vDegreeTmp = new int [vNum];
    memset(vDegreeTmp, 0, sizeof(int) *vNum);
    for (e = 0; e < eNum; e++) {
        v1 = edges[e].v1;
        v2 = edges[e].v2;

        vEdges[v1][vDegreeTmp[v1]] = e;
        vEdges[v2][vDegreeTmp[v2]] = e;

        vAdj[v1][vDegreeTmp[v1]] = v2;
        vAdj[v2][vDegreeTmp[v2]] = v1;

        vDegreeTmp[v1]++;
        vDegreeTmp[v2]++;
    }
    delete[] vDegreeTmp;

    return 1;
}

void prep_remove_from_neighbours(unsigned v) {
    //'remove' v from vAdj: move v to the end part of v_adj's vAdj
    for (unsigned i = 0; i < dVDegree[v]; i++) {
        unsigned adj = vAdj[v][i];
        unsigned indexOfLast = dVDegree[adj] - 1;
        dVDegree[adj] = indexOfLast;
        //vAdj
        for (unsigned j = 0; j <= indexOfLast; j++) {
            if (vAdj[adj][j] == v) {
                vAdj[adj][j] = vAdj[adj][indexOfLast];
                vAdj[adj][indexOfLast] = v;

                unsigned e = vEdges[adj][j];
                vEdges[adj][j] = vEdges[adj][indexOfLast];
                vEdges[adj][indexOfLast] = e;
                break;
            }
        }
    }
}

inline bool prep_handle_degree_0(unsigned v) {
    vInBest[v] = NIN_ABS;
    total_NIN_ABS++;
    return true;
}

inline bool prep_handle_degree_1(unsigned v) {
    // set adj as an necessary vertice
    unsigned adj = vAdj[v][0];
    vInBest[adj] = IN_ABS;
    total_IN_ABS++;

    // !important the sequence cannot be changed
    prep_remove_from_neighbours(adj); //remove adj from all its neighbours

    // set v as an unnecessary vertice
    vInBest[v] = NIN_ABS;
    total_NIN_ABS++;

    //set adjAdj's checkAdjDomV to true
    for (unsigned k = 0, dadj = dVDegree[adj]; k < dadj; k++) {
        unsigned adjAdj = vAdj[adj][k];
        checkAdjDomV[adjAdj] = IN;
    }

    //remove adj's neighbours from adj
    dVDegree[adj] = 0;
    return true;
}

inline bool prep_handle_degree_2(unsigned v) {
    //set u & w as necessary vertexes
    unsigned u = vAdj[v][0];
    unsigned w = vAdj[v][1];
    bool found = false;
    if (dVDegree[u] < dVDegree[w]) {
        for (unsigned i = 0; i < dVDegree[u]; i++) {
            if (vAdj[u][i] == w) {
                found = true;
                break;
            }
        }
    } else {
        for (unsigned i = 0; i < dVDegree[w]; i++) {
            if (vAdj[w][i] == u) {
                found = true;
                break;
            }
        }
    }
    if (found == false) {
        checkAdjDomV[v] = NIN;
        return false;
    }

    vInBest[u] = IN_ABS;
    prep_remove_from_neighbours(u);
    //set adjAdj's checkAdjDomV to true
    for (unsigned k = 0, dadj = dVDegree[u]; k < dadj; k++) {
        unsigned adjAdj = vAdj[u][k];
        checkAdjDomV[adjAdj] = IN;
    }
    dVDegree[u] = 0;//remove u's neighbours

    vInBest[w] = IN_ABS;
    prep_remove_from_neighbours(w);
    //set adjAdj's checkAdjDomV to true
    for (unsigned k = 0, dadj = dVDegree[w]; k < dadj; k++) {
        unsigned adjAdj = vAdj[w][k];
        checkAdjDomV[adjAdj] = IN;
    }
    dVDegree[w] = 0;//remove w's neighbours
    total_IN_ABS = total_IN_ABS + 2;

    //set v as an unnecessary vertex
    vInBest[v] = NIN_ABS;
    total_NIN_ABS++;

    return true;
}

inline bool prep_handle_dom(unsigned v, bool &found) {
    for (unsigned i = 0; i < dVDegree[v]; i++) {
        unsigned adj = vAdj[v][i];
        if (dVDegree[v] > dVDegree[adj])
            continue;

        //check if N[v] <= N[adj]
        bool adjDomV = true;
        for (unsigned j = 0; j < dVDegree[v]; j++) {
            unsigned vadj = vAdj[v][j];
            if (vadj == adj)
                continue;

            bool contain = false;
            for (unsigned k = 0; k < dVDegree[adj]; k++) {
                if (vAdj[adj][k] == vadj) {
                    contain = true;
                    break;
                }
            }
            if (contain == false) {
                adjDomV = false;
                break;
            }
        }

        if (adjDomV) { //v dominates adj or adj dominates v or they dominate each other
            vInBest[adj] = IN_ABS;
            total_IN_ABS++;

            //remove vForSure from its neighbours
            prep_remove_from_neighbours(adj);

            //set adjAdj's checkAdjDomV to true
            for (unsigned k = 0, dadj = dVDegree[adj]; k < dadj; k++) {
                unsigned adjAdj = vAdj[adj][k];
                checkAdjDomV[adjAdj] = IN;
            }

            //remove vForSure's neighbours
            dVDegree[adj] = 0;
            found = true;
            return false;
        }
    }
    checkAdjDomV[v] = NIN;
    return false;
}


void preprocess() {
#ifdef WIN32
    start = clock();
#else
    times(&start);
#endif
    checkAdjDomV = vInC;
    dVDegree = new unsigned[vNum];//static_cast<unsigned *>(static_cast<void *>(dscore));//dynamic vertex degree
    memcpy(dVDegree, vDegree, sizeof (int)*vNum);

    uvSet = new unsigned[vNum];//cVertexes; //untested vertice set
    uvSize = 0; //untested vertice set size

    for (unsigned v = 0; v < vNum; v++) {
        vInBest[v] = UNKNOWN;
        uvSet[v] = v;
        checkAdjDomV[v] = IN;
    }
    uvSize = vNum;

    for (bool found = true; found;) {
        found = false;
        for (unsigned i = 0; i < uvSize; ) {
            unsigned v = uvSet[i];

            bool remove = false;
            if (vInBest[v] != UNKNOWN)
                remove = true;
            else if (dVDegree[v] == 0) {
                remove = prep_handle_degree_0(v);
            } else if (dVDegree[v] == 1) {
                remove = prep_handle_degree_1(v);
            } else if (dVDegree[v] == 2 && checkAdjDomV[v] == IN) {
                remove = prep_handle_degree_2(v);
            }
            if (remove) { //if vForSure in uvSet, remove it
                found = true;
                uvSize--;
                unsigned lastV = uvSet[uvSize];
                uvSet[i] = lastV;
            } else
                i++;
        }
    }

    for (bool findDomHappened = true; findDomHappened;) {
        findDomHappened = false;
        for (unsigned i = 0; i < uvSize; ) {
            unsigned v = uvSet[i];

            bool remove = false;

            if (vInBest[v] != UNKNOWN)
                remove = true;
            else if (dVDegree[v] == 0)
                remove = prep_handle_degree_0(v);
            else if (dVDegree[v] == 1)
                remove = prep_handle_degree_1(v);
            else if (checkAdjDomV[v] == IN)
                remove = prep_handle_dom(v, findDomHappened);

            if (remove) { //if vForSure in uvSet, remove it
                findDomHappened = true;
                uvSize--;
                unsigned lastV = uvSet[uvSize];
                uvSet[i] = lastV;
            } else
                i++;
        }
    }

    //fix uncovEdges
    ueSize = 0;
    for (unsigned i = 0; i < uvSize; i++) {
        unsigned v = uvSet[i];
        for (unsigned j = 0; j < dVDegree[v]; j++) {
            unsigned e = vEdges[v][j];
            unsigned eIndex = uncovEdgeIndexes[e];
            if (eIndex >= ueSize || uncovEdges[eIndex] != e) { //uncovEdges doesn't contains e
                uncovEdges[ueSize] = e;
                uncovEdgeIndexes[e] = ueSize;
                ueSize++;
            }
        }
    }

    //fix cVertices, layout: IN_ABS vertexes | unused area
    cvSize = 0;
    for (unsigned v = 0; v < vNum; v++) {
        if (vInBest[v] == UNKNOWN)
            vInBest[v] = NIN;
        else if (vInBest[v] == IN_ABS) {
            cVertexes[cvSize] = v;
            cVertexIndexes[v] = cvSize;
            cvSize++;
        }
    }

    bestSize = total_IN_ABS;

#ifdef WIN32
    finish = clock();
    double preprocessTime = (static_cast<double>(finish - start)) / CLOCKS_PER_SEC;
#else
    times(&finish);
    double preprocessTime = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime) / sysconf(_SC_CLK_TCK);
#endif
    preprocessTime = round(preprocessTime * 100) / 100.0;

    cout << "c preprocess: total_IN_ABS=" << total_IN_ABS << ", total_NIN_ABS=" << total_NIN_ABS  << " time=" << preprocessTime << endl;

    cout << "c graph: |V'|=" << uvSize << " |E'|=" << ueSize << endl;
}

void free_memory() {
    for (unsigned int v = 0; v < vNum; v++) {
        delete[] vAdj[v];
        delete[] vEdges[v];
    }
    delete[] vInBest;
    delete[] cVertexIndexes;
    delete[] cVertexes;
    delete[] vInC;
    delete[] vDegree;
    delete[] vAdj;
    delete[] vEdges;
    delete[] vTimestamp;
    delete[] dscore;
    delete[] uncovEdgeIndexes;
    delete[] uncovEdges;
    delete[] edges;
    delete[] eTimestamp;
    delete[] uvSet;
    delete[] dVDegree;
}

void reset_remove_cand() {
    unsigned int v, j;
    j = 0;
    for (v = 0; v < vNum; v++) {
        if (vInC[v] == true) {
            cVertexes[j] = v;
            cVertexIndexes[v] = j;
            j++;
        } else cVertexIndexes[v] = 0;
    }

    cvSize = j;
}

void remove_one_vertex_from_c() {
    //TODO to very large graph this function is too time-consuming
    unsigned int best_remove_v = cVertexes[0];
    int    best_dscore = dscore[best_remove_v];

    if (best_dscore != 0) {
        for (unsigned int i = 1; i < cvSize; ++i) {
            unsigned v = cVertexes[i];

            if (dscore[v] == 0) {
                vInBest[v] = NIN;
                bestSize = bestSize - 1;

                remove(v);
                unsigned int last_v_in_c = cVertexes[--cvSize];
                unsigned int index = cVertexIndexes[v];
                cVertexes[index] = last_v_in_c;
                cVertexIndexes[last_v_in_c] = index;
                break;
            } else if (dscore[v] == 1) {
                best_remove_v = v;
                break;
            }

            if (dscore[v] > best_dscore) {
                best_remove_v = v;
                best_dscore = dscore[best_remove_v];
            }
        }
    }

    remove(best_remove_v);

    //remove one vertex from c, and move the last vertex in c to the position
    unsigned int last_v_in_c = cVertexes[--cvSize];
    unsigned int index = cVertexIndexes[best_remove_v];
    cVertexes[index] = last_v_in_c;
    cVertexIndexes[last_v_in_c] = index;
}


//update the best vertex in C
static unsigned int cand_count = 50;
unsigned int choose_remove_v() {
    unsigned int i, v;

    unsigned int best_v = cVertexes[static_cast<unsigned int>(pcg32_fast()) % cvSize];
    if (static_cast<unsigned int>(pcg32_fast()) % 10 < 9) {
        for (i = 1; i < cand_count; ++i) {
            v = cVertexes[static_cast<unsigned int>(pcg32_fast()) % cvSize];

            if (dscore[v] < dscore[best_v])
                continue;
            else if (dscore[v] > dscore[best_v])
                best_v = v;
            else if (vTimestamp[v] < vTimestamp[best_v]) {
                best_v = v;
            }
        }
    }

    return best_v;
}



inline
void uncover(unsigned int e) {
    uncovEdgeIndexes[e] = ueSize;
    uncovEdges[ueSize++] = e;
    eTimestamp[e] = step;
}


inline
void cover(unsigned int e) {
    unsigned int index, last_uncov_edge;

    //since the edge is satisfied, its position can be reused to store the last_uncov_edge
    last_uncov_edge = uncovEdges[--ueSize];
    index = uncovEdgeIndexes[e];
    uncovEdges[index] = last_uncov_edge;
    uncovEdgeIndexes[last_uncov_edge] = index;
}

void init_sol() {
    //use dVDegree as vDegree
    int *dscoreTemp = new int[vNum];
    VerticeStatus *vInCTemp = new VerticeStatus[vNum];
    unsigned int cvSizeTemp = 0;

    memset(vTimestamp, 0, sizeof(unsigned int) * vNum);

    init_vertex_greedy_vc(dscoreTemp, vInCTemp, cvSizeTemp);
    check_solution(vInCTemp, cvSizeTemp + total_IN_ABS, dscoreTemp);
    std::swap(vInC, vInCTemp);
    std::swap(dscore, dscoreTemp);
    cvSize = cvSizeTemp;

    init_edge_greedy_vc(dscoreTemp, vInCTemp, cvSizeTemp);
    check_solution(vInCTemp, cvSizeTemp + total_IN_ABS, dscoreTemp);
    if (cvSizeTemp < cvSize) {
        std::swap(vInC, vInCTemp);
        std::swap(dscore, dscoreTemp);
        cvSize = cvSizeTemp;
    }

    init_match_vc(dscoreTemp, vInCTemp, cvSizeTemp);
    check_solution(vInCTemp, cvSizeTemp + total_IN_ABS, dscoreTemp);
    if (cvSizeTemp < cvSize) {
        std::swap(vInC, vInCTemp);
        std::swap(dscore, dscoreTemp);
        cvSize = cvSizeTemp;
    }

    //construct other structures
    update_best_sol();

    //cAllVertexes layout:IN_ABS vertices | IN vertices | unused area
    //cVertexes layout:IN vertices | unused area
    cAllVertexes = cVertexes;
    cVertexes = &(cAllVertexes[total_IN_ABS]);
    cAVSize = cvSize + total_IN_ABS;
    cvSize = 0;
    for (unsigned i = 0; i < uvSize; i++) {
        unsigned v = uvSet[i];
        if (in(vInC, v)) {
            cVertexes[cvSize] = v;
            cVertexIndexes[v] = cvSize;
            cvSize++;
        }
    }

    ueSize = 0;

    cout << "c init best = " << cAVSize << endl;
    delete [] dscoreTemp;
    delete [] vInCTemp;
}

void init_vertex_greedy_vc(int *dscore1, VerticeStatus *vInC1, unsigned int &cvSize1) {  //super fast
#ifdef WIN32
    start = clock();
#else
    times(&start);
#endif
    struct ScoreVertex {
        int score;
        unsigned int vertex;
    };

    function<bool(const ScoreVertex &, const ScoreVertex &)> compFunc = [](const ScoreVertex & v1, const ScoreVertex & v2) {
        return v1.score < v2.score;//<:大的在前面，小的在后面
    };

    ScoreVertex *candidates = new ScoreVertex[uvSize];
    unsigned int candSize = 0;
    cvSize1 = 0;
    unsigned int ueSize1 = ueSize;

    unsigned *readdCandVertices = new unsigned[uvSize];
    unsigned racSize = 0;
    unsigned *highScoreVertices = new unsigned[uvSize];
    unsigned hsvSize = 0;
    int highestScore = static_cast<int>(vNum);

    for (unsigned int v = 0; v < vNum; v++) {
        vInC1[v] = vInBest[v];
        dscore1[v] = static_cast<int>(dVDegree[v]);
        if (vInC1[v] != IN_ABS && vInC1[v] != NIN_ABS) {
            vInC1[v] = NIN;
            candidates[candSize].vertex = v;
            candidates[candSize].score = dscore1[v];
            candSize++;
        }
    }

    make_heap(candidates, candidates + uvSize, compFunc);

    for (; candSize > 0 && ueSize1 > 0;) {
        highestScore = candidates[0].score;
        int gainValue;
        unsigned int v;
        //pop all vertices with bigest score value
        while (candSize > 0) {
            gainValue = candidates[0].score;
            v = candidates[0].vertex;
            if (gainValue == highestScore) {
                pop_heap(candidates, candidates + candSize, compFunc);
                candSize--;
                if (gainValue != dscore1[v]) {//readd the vertexes whose gain value has changed
                    if (dscore1[v] > 0) {
                        readdCandVertices[racSize] = v;
                        racSize++;
                    }
                } else {//add to highScoreVertices
                    highScoreVertices[hsvSize] = v;
                    hsvSize++;
                }
            } else {
                break;
            }
        }

        //Randomly select and add
        while (hsvSize > 0) {
            unsigned randomIndex = pcg32_fast() % hsvSize;
            hsvSize--;
            v = highScoreVertices[randomIndex];
            highScoreVertices[randomIndex] = highScoreVertices[hsvSize];

            if (dscore1[v] != highestScore) {
                if (dscore1[v] > 0) {
                    readdCandVertices[racSize] = v;
                    racSize++;
                }
            } else {// Add v to c
                cvSize1++;
                dscore1[v] = dscore1[v] * -1;
                vInC1[v] = IN;

                unsigned int edgeCount = dVDegree[v];
                for (unsigned int i = 0; i < edgeCount; ++i) {
                    unsigned int n = vAdj[v][i];//v's i'th neighbor

                    if (not_in(vInC1, n)) { //this adj isn't in cover set
                        dscore1[n]--;
                        ueSize1--;
                    } else {
                        dscore1[n]++;
                    }
                }
            }
        }

        //readd
        while (racSize > 0) {
            racSize--;
            v = readdCandVertices[racSize];

            candidates[candSize].vertex = v;
            candidates[candSize].score = dscore1[v];
            push_heap(candidates, candidates + candSize, compFunc);
            candSize++;
        }
    }

    for (unsigned i = 0; i < uvSize; i++) {
        unsigned v = uvSet[i];
        if (dscore1[v] == 0 && in(vInC1, v)) {
            cvSize1--;
            // remove v from c
            vInC1[v] = NIN;

            unsigned int edgeCount = dVDegree[v];
            for (unsigned int i = 0; i < edgeCount; ++i) {
                unsigned int n = vAdj[v][i];//v's i'th neighbor
                if (not_in(vInC1, n))
                    dscore1[n]++;
                else
                    dscore1[n]--;
            }
        }
    }
    delete [] readdCandVertices;
    delete [] highScoreVertices;
    delete [] candidates;

#ifdef WIN32
    finish = clock();
    double initTime = (static_cast<double>(finish - start)) / CLOCKS_PER_SEC;
#else
    times(&finish);
    double initTime = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime) / sysconf(_SC_CLK_TCK);
#endif

    initTime = round(initTime * 100) / 100.0;

    cout << "c init vertex_greedy_vc=" << cvSize1 << "/" << uvSize << "/" << (cvSize1 + total_IN_ABS) << " time=" << initTime << endl;
}

void init_edge_greedy_vc(int *dscore2, VerticeStatus *vInC2, unsigned int &cvSize2) {
#ifdef WIN32
    start = clock();
#else
    times(&start);
#endif

    unsigned int v, e;
    unsigned int v1, v2;

    memcpy(vInC2, vInBest, vNum);
    memset(dscore2, 0, sizeof(int) * vNum);
    cvSize2 = 0;

    for (unsigned i = 0; i < ueSize; i++) {
        e = uncovEdges[i];
        v1 = edges[e].v1;
        v2 = edges[e].v2;

        if (not_in(vInC2, v1) && not_in(vInC2, v2)) { //if uncovered, choose the endpoint with higher degree
            if (dVDegree[v1] > dVDegree[v2]) {
                vInC2[v1] = IN;
            } else {
                vInC2[v2] = IN;
            }
            cvSize2++;
        }
    }

    //calculate dscores: for the vertexes in C, reduce its score by 1; for the vertexes not in C, remain its score to 0, including IN_ABS and NIN_ABS vertexes
    for (unsigned i = 0; i < ueSize; i++) {
        e = uncovEdges[i];
        v1 = edges[e].v1;
        v2 = edges[e].v2;

        if (in(vInC2, v1) && not_in(vInC2, v2)) dscore2[v1]--;
        else if (in(vInC2, v2) && not_in(vInC2, v1)) dscore2[v2]--;
    }

    //remove redundent vertices
    for (unsigned i = 0; i < uvSize; i++) {
        v = uvSet[i];
        if (in(vInC2, v) && dscore2[v] == 0) {
            vInC2[v] = NIN;
            dscore2[v] = -1 * dscore2[v];
            unsigned int degree = dVDegree[v];
            for (unsigned int i = 0; i < degree; i++) {
                unsigned int vadj = vAdj[v][i];
                if (vInC2[vadj] == false) {
                    dscore2[vadj]++;
                } else
                    dscore2[vadj]--;
            }
            cvSize2--;
        }
    }

#ifdef WIN32
    finish = clock();
    double initTime = (static_cast<double>(finish - start)) / CLOCKS_PER_SEC;
#else
    times(&finish);
    double initTime = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime) / sysconf(_SC_CLK_TCK);
#endif
    initTime = round(initTime * 100) / 100.0;

    cout << "c init edge_greedy_vc=" << cvSize2 << "/" << uvSize << "/" << (cvSize2 + total_IN_ABS) << " time=" << initTime << endl;
}

void init_match_vc(int *dscore3, VerticeStatus *vInC3, unsigned int &cvSize3) {
#ifdef WIN32
    start = clock();
#else
    times(&start);
#endif

    unsigned int v, e;
    unsigned int v1, v2;

    memset(dscore3, 0, sizeof(int) * vNum);
    memcpy(vInC3, vInBest, vNum);

    cvSize3 = 0;
    for (unsigned i = 0; i < ueSize; i++) {
        e = uncovEdges[i];
        v1 = edges[e].v1;
        v2 = edges[e].v2;

        if (not_in(vInC3, v1) && not_in(vInC3, v2)) {
            vInC3[v1] = IN;
            vInC3[v2] = IN;
            cvSize3 = cvSize3 + 2;
        }
    }

    //calculate dscores
    for (unsigned i = 0; i < ueSize; i++) {
        e = uncovEdges[i];
        v1 = edges[e].v1;
        v2 = edges[e].v2;

        if (in(vInC3, v1) && not_in(vInC3, v2))
            dscore3[v1]--;
        else if (in(vInC3, v2) && not_in(vInC3, v1))
            dscore3[v2]--;
    }

    //remove redundent vertices
    for (unsigned i = 0; i < uvSize; i++) {
        v = uvSet[i];
        if (in(vInC3, v) && dscore3[v] == 0) {
            vInC3[v] = NIN;
            dscore3[v] = -1 * dscore3[v];
            unsigned int degree = dVDegree[v];
            for (unsigned int i = 0; i < degree; i++) {
                unsigned int vadj = vAdj[v][i];
                if (not_in(vInC3, vadj)) {
                    dscore3[vadj]++;
                } else
                    dscore3[vadj]--;
            }
            cvSize3--;
        }
    }

#ifdef WIN32
    finish = clock();
    double initTime = (static_cast<double>(finish - start)) / CLOCKS_PER_SEC;
#else
    times(&finish);
    double initTime = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime) / sysconf(_SC_CLK_TCK);
#endif
    initTime = round(initTime * 100) / 100.0;

    cout << "c init match_vc=" << cvSize3 << "/" << uvSize << "/" << (cvSize3 + total_IN_ABS) << " time=" << initTime << endl;
}


void add(unsigned int v) {
    vInC[v] = IN;
    dscore[v] = -dscore[v];

    unsigned int i, e, n;

    unsigned int edgeCount = dVDegree[v];

    for (i = 0; i < edgeCount; ++i) {
        e = vEdges[v][i];// v's i'th edge
        n = vAdj[v][i];//v's i'th neighbor

        if (not_in(vInC, n)) { //this adj isn't in cover set
            dscore[n]--;
            //conf_change[n] = 1;

            cover(e);
        } else {
            dscore[n]++;
        }
    }

}

void remove(unsigned int v) {
    vInC[v] = NIN;
    dscore[v] = -dscore[v];

    unsigned int i, e, n;

    unsigned int edge_count = dVDegree[v];
    for (i = 0; i < edge_count; ++i) {
        e = vEdges[v][i];
        n = vAdj[v][i];

        if (not_in(vInC, n)) { //this adj isn't in cover set
            dscore[n]++;
            uncover(e);
        } else {
            dscore[n]--;
        }
    }
}

/*On solution*/

void print_solution() {
    for (unsigned int i = 0; i < vNum; i++) {
        if (in(vInBest, i)) //output vertex cover
            cout << (i + 1) << '\t';
    }
    cout << endl;
}


int check_solution() {
    for (unsigned int e = 0; e < eNum; ++e) {
        if (!in(vInBest, edges[e].v1) && !in(vInBest, edges[e].v2)) {
            cout << "c error: uncovered edge " << e << endl;
            return 0;
        }
    }

    unsigned int verified_vc_size = 0;
    for (unsigned int i = 0; i < vNum; i++) {
        if (in(vInBest, i))
            verified_vc_size++;
    }

    if (bestSize == verified_vc_size)
        return 1;
    else {
        cout << "c error: claimed best found vc size!=verified vc size" << endl;
        cout << "c claimed best found vc size=" << bestSize << endl;
        cout << "c verified vc size=" << verified_vc_size << endl;
        return 0;
    }
}

bool check_solution(VerticeStatus *vInCTest, unsigned int cvSizeTest, int *scoreTest) {
    bool passAll = true;
    //test if cover all edges
    for (unsigned int e = 0, c = 0; e < eNum && c < 5; ++e) {
        if (not_in(vInCTest, edges[e].v1) && not_in(vInCTest, edges[e].v2)) {
            cout << "c error: uncovered edge " << e << endl;
            passAll = false;
            c++;
        }
    }

    //test if the claimed size is right
    unsigned int verified_vc_size = 0;
    for (unsigned int i = 0; i < vNum; i++) {
        if (in(vInCTest, i))
            verified_vc_size++;
    }
    if (cvSizeTest != verified_vc_size) {
        cout << "c error: claimed vc size(" << cvSizeTest << ") != " <<
             "verified vc size(" << verified_vc_size << ")" << endl;
        passAll = false;
    }

    //test if score is right
    if (scoreTest != nullptr) {
        int *scoreReal = new int[vNum];
        //calculate dscores
        for (unsigned int v = 0; v < vNum; v++) {
            if (vInCTest[v] == IN_ABS || vInCTest[v] == NIN_ABS)
                scoreReal[v] = 0;
            else {
                int factor = in(vInCTest, v) ? -1 : 1;

                scoreReal[v] = 0;
                unsigned int edgeCount = vDegree[v];
                for (unsigned int i = 0; i < edgeCount; i++) {
                    unsigned int adj = vAdj[v][i];
                    if (not_in(vInCTest, adj))
                        scoreReal[v]++;
                }
                scoreReal[v] = scoreReal[v] * factor;
            }
        }

        //compare
        for (unsigned int v = 0; v < vNum; v++) {
            if (scoreReal[v] != scoreTest[v]) {
                cout << "c error: for vertex " << v << ": real score(" << scoreReal[v] << ") != claimed score(" << scoreTest[v] << ")" << "   degree=" << vDegree[v] << endl;
                passAll = false;
                break;
            }
        }

        delete [] scoreReal;
    }

    return passAll;
}

inline bool in(VerticeStatus *vs, unsigned v) {
    if (vs[v] > UNKNOWN) { // vs[v] == IN || vs[v] == IN_ABS
        return true;
    } else {
        return false;
    }
}

inline bool not_in(VerticeStatus *vs, unsigned v) {
    if (vs[v] < UNKNOWN) { //vs[v] == NIN || vs[v] == NIN_ABS
        return true;
    } else {
        return false;
    }
}

#endif // EAVC_H
