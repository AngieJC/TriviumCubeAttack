#include <vector>
#include <map>
#include <algorithm>
#include <set>
#include <iostream>
#include <bitset>
#include <time.h>
#include <cmath>
#include <chrono>
#include <malloc.h>
#include <pthread.h>

#include "gurobi_c++.h"

#define TimeLimit 0         // 0 for no TimeLimit; 1 for dynamic TimeLimit; others for const TimeLimit
#define TimeLimitBase 5
#define PRINT_GRAPH 0
#define SAVE_MILP 0

#define PATTERN2 1

#define NEWPATTERN 0

#define PATTERNA 0
#define PATTERNB 0
#define PATTERNC 0
#define PATTERND 0

#define PATTERNE 0
#define PATTERNF 0

#define PATTERNG 0

#define PATTERNH 0
#define PATTERNI 0

#define PATTERNJ 0
#define PRINTEQU 1

#define DEVIDEUPBOUND 11
#define MULTITHREAD 32

using namespace std;
using namespace std::chrono;

void printEq(vector<set<int>> const& v) {
    if (v.empty()) cout << " 0 " << endl;
    else {
        for (auto const& s : v) {
            if (s.empty()) cout << " 1 ";
            else {
                for (auto x : s) cout << "v" << x << ".";
            }
            cout << " + ";
        }
        cout << endl;
    }
}

struct cmpBitset288 {
    bool operator()(const bitset<288>& a, const bitset<288>& b) const {
        for (int i = 0; i < 288; i++) {
            if (a[i] < b[i])
                return true;
            else if (a[i] > b[i])
                return false;
        }
        return false;
    }
};

int display_result_anf(set<int> cube, map<bitset<288>, int, cmpBitset288>& countingBox) {
    int i;
    auto it2 = countingBox.begin();
    while (it2 != countingBox.end()) {
        if (((*it2).second % 2) == 1) {
            cout << ((*it2).second % 2) << " | " << (*it2).second << "\t";
            bitset<288> tmp = (*it2).first;
            for (i = 0; i < 80; i++) if ((tmp[i] == 1))
                cout << "k" << (i) << " ";
            /*for (i = 0; i < 80; i++) if ((tmp[93 + i] == 1))
              cout << "v" << (i) << " ";
            // Uncomment for the cube in the superpoly.*/
            for (i = 0; i < 80; i++) if ((tmp[93 + i] == 1 && cube.find(i) != cube.end()))//!binary_search(cube.begin(),cube.end(),i))) 
                cout << "v" << (i) << " ";
            cout << endl;
        }
        it2++;
    }
    return 0;
}


void printSys(map<int, vector<set<int>>> const& sys) {
    for (auto const& p : sys) {
        cout << "v" << p.first << " = ";
        printEq(p.second);
    }
}


int getVarTrivium(int x) {
    if (x < 288) return x;
    auto xmod = x % 288;
    if (xmod == 0 || xmod == 93 || xmod == 177) return x;
    return getVarTrivium(x - 288 - 1);
}

map<int, vector<set<int>>> generateTrivium(int R) {
    map<int, vector<set<int>>> sys;
    int root = 288 * (R + 1);
    auto& root_eq = sys[root];
    vector<int> indexes_root = { 65, 92, 161, 176, 242, 287 };
    set<int> waiting;
    for (auto i : indexes_root) {
        root_eq.emplace_back(set<int>({ getVarTrivium(288 * R + i) }));
        waiting.emplace(getVarTrivium(288 * R + i));
    }
    while (!waiting.empty()) {
        auto x = *(waiting.rbegin());
        waiting.erase(x);
        if (x / 288 == 0) continue;
        auto& eq = sys[x];
        switch (x % 288) {
        case 0:
            x -= 288;
            eq.emplace_back(set<int>({ getVarTrivium(x + 68) }));
            eq.emplace_back(set<int>({ getVarTrivium(x + 242) }));
            eq.emplace_back(set<int>({ getVarTrivium(x + 287) }));
            eq.emplace_back(set<int>({ getVarTrivium(x + 285), getVarTrivium(x + 286) }));
            waiting.emplace(getVarTrivium(x + 242));
            waiting.emplace(getVarTrivium(x + 287));
            waiting.emplace(getVarTrivium(x + 68));
            waiting.emplace(getVarTrivium(x + 285));
            waiting.emplace(getVarTrivium(x + 286));
            break;
        case 93:
            x -= 288 + 93;
            eq.emplace_back(set<int>({ getVarTrivium(x + 170) }));
            eq.emplace_back(set<int>({ getVarTrivium(x + 65) }));
            eq.emplace_back(set<int>({ getVarTrivium(x + 92) }));
            eq.emplace_back(set<int>({ getVarTrivium(x + 90), getVarTrivium(x + 91) }));
            waiting.emplace(getVarTrivium(x + 65));
            waiting.emplace(getVarTrivium(x + 92));
            waiting.emplace(getVarTrivium(x + 170));
            waiting.emplace(getVarTrivium(x + 90));
            waiting.emplace(getVarTrivium(x + 91));
            break;
        case 177:
            x -= 288 + 177;
            eq.emplace_back(set<int>({ getVarTrivium(x + 263) }));
            eq.emplace_back(set<int>({ getVarTrivium(x + 161) }));
            eq.emplace_back(set<int>({ getVarTrivium(x + 176) }));
            eq.emplace_back(set<int>({ getVarTrivium(x + 174), getVarTrivium(x + 175) }));
            waiting.emplace(getVarTrivium(x + 161));
            waiting.emplace(getVarTrivium(x + 176));
            waiting.emplace(getVarTrivium(x + 263));
            waiting.emplace(getVarTrivium(x + 174));
            waiting.emplace(getVarTrivium(x + 175));
            break;
        }
    }
    return sys;
}

map<int, vector<int>> reverseSys(map<int, vector<set<int>>> const& sys) {
    map<int, vector<int>> res;
    for (auto const& eq : sys) {
        for (auto const& s : eq.second) {
            for (auto x : s) res[x].emplace_back(eq.first);
        }
    }
    return res;
}

int getT(int x, map<int, vector<set<int>>> const& sys, map<int, int>& T, map<int, int>& TT);
int getT(set<int> const& s, map<int, vector<set<int>>> const& sys, map<int, int>& T, map<int, int>& TT);

int getTT(int x, map<int, vector<set<int>>> const& sys, map<int, int>& T, map<int, int>& TT) {
    auto it = TT.find(x);
    if (it != TT.end()) return it->second;
    auto it_eq1 = sys.find(x);
    auto it_eq2 = sys.find(x + 288);
    auto& TTx = TT[x];
    if (it_eq1 == sys.end() && it_eq2 == sys.end()) {
        TTx = getT(x, sys, T, TT) + getT(x + 288, sys, T, TT);
    }
    else {
        auto const& eq1 = (it_eq1 != sys.end()) ? it_eq1->second : vector<set<int>>(1, set<int>({ x }));
        auto const& eq2 = (it_eq2 != sys.end()) ? it_eq2->second : vector<set<int>>(1, set<int>({ x + 288 }));
        TTx = 0;
        for (auto const& s1 : eq1) {
            for (auto s2 : eq2) {
                s2.insert(s1.begin(), s1.end());
                int tmp = getT(s2, sys, T, TT);
                if (tmp > TTx) TTx = tmp;
            }
        }
    }
    return TTx;
}

int getT(set<int> const& s, map<int, vector<set<int>>> const& sys, map<int, int>& T, map<int, int>& TT) {
    if (s.size() == 0) return 0;
    if (s.size() == 1) return getT(*s.begin(), sys, T, TT);
    auto ss = s;
    auto x = *ss.begin();
    ss.erase(ss.begin());
    if (ss.count(x + 288) == 0) return getT(x, sys, T, TT) + getT(ss, sys, T, TT);
    else {
        ss.erase(x + 288);
        return getTT(x, sys, T, TT) + getT(ss, sys, T, TT);
    }
}

int getT(int x, map<int, vector<set<int>>> const& sys, map<int, int>& T, map<int, int>& TT) {
    auto it = T.find(x);
    if (it != T.end()) return it->second;
    auto& Tx = T[x];
    Tx = 0;
    auto it_eq = sys.find(x);
    if (it_eq != sys.end()) {
        for (auto const& s : it_eq->second) {
            int local_best = getT(s, sys, T, TT);
            if (local_best > Tx) Tx = local_best;
        }
    }
    return Tx;
}

void triviumCore(GRBModel& model, vector<GRBVar>& x, int i1, int i5, int i2, int i3, int i4)
{
    int Ineq[][11] = {
  {0, -1, -1, 0, -1, -1, 1, 1, 0, 1, 1},
  {0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0},
  {0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0},
  {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1},
  {0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0},
  {0, -1, 0, -1, -1, -1, 1, 0, 1, 1, 1},
  {0, 0, -1, 1, 0, 0, 0, 1, 0, 0, 0},
  {0, 0, 1, -1, 0, 0, 0, 0, 1, 0, 0},
  {2, 0, 1, 0, 1, 0, -1, 0, 0, -1, -1},
  {0, 1, 1, 0, 1, 1, 0, 0, 0, -1, 0},
  {3, 1, 0, 0, 1, 0, 0, -1, -1, -1, -1},
  {2, 0, 0, 1, 1, 0, -1, 0, 0, -1, -1},
  {0, 1, 0, 1, 1, 1, 0, 0, 0, -1, 0},
  {3, 0, 0, 0, 1, 1, -1, -1, -1, -1, 0}
    };
    GRBVar y1 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar y2 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar y3 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar y4 = model.addVar(0, 1, 0, GRB_BINARY);
    GRBVar y5 = model.addVar(0, 1, 0, GRB_BINARY);

    for (auto it : Ineq)
        model.addConstr(it[0] + it[1] * x[i1] + it[2] * x[i2] + it[3] * x[i3] + it[4] * x[i4] + it[5] * x[i5] + it[6] * y1 + it[7] * y2 + it[8] * y3 + it[9] * y4 + it[10] * y5 >= 0);

    x[i1] = y1;
    x[i2] = y2;
    x[i3] = y3;
    x[i4] = y4;
    x[i5] = y5;
}

int BackExpandPolynomial(int rounds, vector<bitset<288> >& term)
{
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_PoolSearchMode, 2);
    env.set(GRB_IntParam_PoolSolutions, 2000000000);

    GRBModel model = GRBModel(env);

    vector<GRBVar> s(288);
    for (int i = 0; i < 288; i++)
        s[i] = model.addVar(0, 1, 0, GRB_BINARY);

    vector<GRBVar> works = s;
    for (int r = 0; r < rounds; r++)
    {
        triviumCore(model, works, 65, 170, 90, 91, 92);
        triviumCore(model, works, 161, 263, 174, 175, 176);
        triviumCore(model, works, 242, 68, 285, 286, 287);

        vector<GRBVar> temp = works;
        for (int i = 0; i < 288; i++)
            works[(i + 1) % 288] = temp[i];
    }

    GRBLinExpr nk = 0;
    for (int i = 0; i < 288; i++)
        if ((i == 65) || (i == 92) || (i == 161) || (i == 176) || (i == 242) || (i == 287))
            nk += works[i];
        else
            model.addConstr(works[i] == 0);
    model.addConstr(nk == 1);

    model.update();
    model.optimize();

    map<bitset<288>, int, cmpBitset288> counterMap;

    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
    {
        double time = model.get(GRB_DoubleAttr_Runtime);
        cout << "Time Used: " << time << "sec" << endl;

        int solCount = model.get(GRB_IntAttr_SolCount);
        cout << "Raw Solutions: " << solCount << endl;

        bitset<288> start;
        for (int i = 0; i < solCount; i++)
        {
            model.set(GRB_IntParam_SolutionNumber, i);

            for (int j = 0; j < 288; j++)
                if (round(s[j].get(GRB_DoubleAttr_Xn)) == 1)
                    start[j] = 1;
                else
                    start[j] = 0;
            counterMap[start]++;
        }
    }
    else if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
    {
        cout << "No terms " << endl;
        exit(-2);
    }
    else
    {
        cout << "Other status " << GRB_IntAttr_Status << endl;
        exit(-1);
    }

    for (auto it : counterMap)
        if (it.second % 2 == 1)
            term.push_back(it.first);
    cout << "Exact Solutions: " << term.size() << endl;
    return 0;
}

void bitset2roots(vector<bitset<288>>& term, vector<set<int>>& vsetroot, int zRound, int backRound) {
    int round = zRound - backRound, index = 0;
    for (auto monomial : term) {
        set<int> temp;
        for (int i = 0; i < 288; i++) {
            if (monomial[i]) {
                if (i < 93) {
                    index = 0;
                }
                else if (i < 177) {
                    index = 93;
                }
                else {
                    index = 177;
                }
                temp.emplace((round - i + index) * 288 + index);
            }
        }
        vsetroot.emplace_back(temp);
    }
}

// multi threads
struct args4MultiThread {
    int UID, rootSize, vcubesize;
    int* setrootarray;
    // set<int>** subnodearray;
    vector<set<int>>* vsetroot;
    map<int, int>* T;
    map<int, int>* TT;
    map<int, vector<set<int>>>* sys;
    map<int, vector<set<int>>>* sys2;
    pthread_mutex_t* mutex;
    pthread_barrier_t* barrier;
};
void* devideMILP(void* ptr) {
	args4MultiThread* args = (args4MultiThread*)(ptr);
	int burpRange = (int)(pow(4, args->rootSize) / MULTITHREAD);
    set<int>** subnodearray = (set<int>**)malloc(4 * sizeof(set<int>*));
    vector<set<int>> newroots;
	for (int i = 0; i < burpRange; i++) {
		int range = i * MULTITHREAD + args->UID;
		set<int> subRoot;
		int highnum = 0;
		for (int j = 0; j < args->rootSize; j++) {
			int k = 0;
			for (auto& subnode : args->sys2->at(args->setrootarray[j])) {
				subnodearray[k++] = &subnode;
			}
			int dim = (int)((range - highnum) / pow(4, args->rootSize - j - 1));
			if (dim < 3) {
				subRoot.insert(*(subnodearray[dim])->begin());
			}
			else {
				subRoot.insert(*(subnodearray[dim])->begin());
				subRoot.insert(*(subnodearray[dim])->rbegin());
			}
			highnum += (int)(dim * pow(4, args->rootSize - j - 1));
		}
		if (getT(subRoot, *args->sys, *args->T, *args->TT) >= args->vcubesize) {
            newroots.emplace_back(subRoot);
		}
	}
    pthread_mutex_lock(args->mutex);
    args->vsetroot->insert(args->vsetroot->end(), newroots.begin(), newroots.end());
    pthread_mutex_unlock(args->mutex);
    vector<set<int>>().swap(newroots);
    malloc_trim(0);
	return NULL;
}

void milpTrivium(int R, set<int> vcube, set<int> v0, int backRound) {
    auto sys = generateTrivium(R);
    auto rsys = reverseSys(sys);

    try {

        // Create an environment
        GRBEnv env = GRBEnv(true);
        // env.set("LogFile", "mip1.log");
        env.set(GRB_IntParam_LogToConsole, 0);
        env.set(GRB_DoubleParam_PoolGap, GRB_INFINITY);
        env.start();

        int root = sys.rbegin()->first;

        {
            set<int> tmp;
            for (auto x : vcube) tmp.emplace(x + 93);
            vcube = move(tmp);
        }

        {
            set<int> tmp;
            for (auto x : v0) tmp.emplace(x + 93);
            v0 = move(tmp);
        }

        //add 0s for state 0 Trivium
        for (int i = 80; i < 93; ++i) v0.emplace(i);
        for (int i = 93 + 80; i < 285; ++i) v0.emplace(i);

        //propagate 0s  (not really needed)
        for (auto const& p : sys) {
            bool flag0 = true;
            for (unsigned i = 0; i < 3; ++i) {
                if (v0.count(*p.second[i].begin()) == 0) flag0 = false;
            }
            auto y1 = *p.second[3].begin();
            auto y2 = *p.second[3].rbegin();
            if (v0.count(y1) == 0 && v0.count(y2) == 0) flag0 = false;
            if (flag0) v0.emplace(p.first);
        }

        //setup dynamic programming for "cost"
        map<int, int> T, TT;
        for (auto y : v0) T[y] = -100;
        for (auto y : vcube) T[y] = 1;

        vector<set<int>> vsetroot;
        vector<set<int>> allvsetroot;
        vector<bitset<288>> midTermsBitset;
        BackExpandPolynomial(backRound, midTermsBitset);
        bitset2roots(midTermsBitset, allvsetroot, R, backRound);
        for (auto x : allvsetroot) {
            if (getT(x, sys, T, TT) >= vcube.size()) {
                vsetroot.emplace_back(x);
            }
        }

        int cptt = 0;
        int nSolutions = 0;
        int deletedModels = 0;
        map<bitset<288>, int, cmpBitset288> countingBox;


        // handle each model
        {
            uint64_t iroot = 0;
            while (iroot < vsetroot.size()) {
                auto setroot = move(vsetroot[iroot++]);

                // generate submodel
                auto waiting = setroot;
                map<int, vector<set<int>>> sys2;
                while (!waiting.empty()) {
                    auto x = *waiting.rbegin();
                    waiting.erase(x);
                    auto it = sys.find(x);
                    if (it == sys.end()) continue;
                    sys2[x] = it->second;
                    for (auto const& s : sys2[x]) waiting.insert(s.begin(), s.end());
                }
                if (sys2.find(288) == sys2.end())
                {
                    continue;
                }

                int Nvar = 0;

                map<pair<int, int>, int> mapVars;
                for (auto const& p : sys2) {
                    for (auto const& s : p.second) {
                        for (auto x : s) mapVars[make_pair(p.first, x)] = Nvar;
                        Nvar++;
                    }
                }

                auto rsys2 = reverseSys(sys2);

                // Create an empty model
                GRBModel* model = new GRBModel(env);

                // Create variables and setup strategy
                vector<GRBVar> X(Nvar);
                for (auto const& p : sys2) {
                    for (auto const& s : p.second) {
                        auto x = *s.begin();
                        {
                            if (s.size() == 1 || getT(s, sys, T, TT) <= 1) X[mapVars[make_pair(p.first, x)]] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY);
                            else X[mapVars[make_pair(p.first, x)]] = model->addVar(0.0, 1.0, 1.0, GRB_BINARY);
                            if (getT(s, sys, T, TT) <= 0) X[mapVars[make_pair(p.first, x)]].set(GRB_IntAttr_BranchPriority, -300);
                            else if (s.size() == 2) X[mapVars[make_pair(p.first, x)]].set(GRB_IntAttr_BranchPriority, getT(s, sys, T, TT));
                        }
                    }
                    for (int i = 0; i < 3; ++i) {
                        X[mapVars[make_pair(p.first, *p.second[i].begin())]].set(GRB_DoubleAttr_VarHintVal, 0.0);
                    }
                }

                // remove some impossible cases using dynamic programming
                {
                    map<int, vector<int>> waiting;
                    waiting[*setroot.begin()] = vector<int>(4, vcube.size());
                    while (!waiting.empty()) {
                        auto pw = *waiting.rbegin();
                        waiting.erase(pw.first);
                        if (pw.first < 288) continue;
                        auto const& eq = sys2.at(pw.first);
                        for (int i = 0; i <= 3; ++i) {
                            auto const& s = eq[i];
                            if (getT(s, sys, T, TT) < pw.second[i]) {
                                //if (getT(s, sys, T, TT) < vcube.size()) {
                                    // impossible to have pw.first --> eq[i]
                                    // model.addConstr(X[mapVars[make_pair(pw.first, *s.begin())]] == 0);
                            }
                            else {
                                if (i < 3) {
                                    auto it = waiting.find(*s.begin());
                                    if (it == waiting.end()) waiting.emplace(*s.begin(), vector<int>(4, pw.second[i]));
                                    else {
                                        for (int j = 0; j < 4; ++j) it->second[j] = min(it->second[j], pw.second[i]);
                                    }
                                }
                                else {
                                    auto x1 = *s.begin();
                                    auto x2 = *s.rbegin();
                                    if (x1 < 288 || x2 < 288) continue;
                                    auto const& eq1 = sys2.at(x1);
                                    auto const& eq2 = sys2.at(x2);
                                    auto it1 = waiting.find(x1);
                                    if (it1 == waiting.end()) it1 = waiting.emplace(x1, vector<int>(4, vcube.size())).first;
                                    auto it2 = waiting.find(x2);
                                    if (it2 == waiting.end()) it2 = waiting.emplace(x2, vector<int>(4, vcube.size())).first;
                                    int m1 = 0;
                                    for (int j = 0; j < 3; ++j) m1 = max(m1, getT(eq1[j], sys, T, TT));
                                    int m2 = 0;
                                    for (int j = 0; j < 3; ++j) m2 = max(m2, getT(eq2[j], sys, T, TT));
                                    if (m1 + m2 < pw.second[i]) {
                                        // At least one of the two branches has to double
                                        model->addConstr(X[mapVars[make_pair(pw.first, *s.begin())]] <= X[mapVars[make_pair(x1, *eq1[3].begin())]] + X[mapVars[make_pair(x2, *eq2[3].begin())]]);
                                    }
                                    for (int j1 = 0; j1 < 4; ++j1) {
                                        for (int j2 = 0; j2 < 4; ++j2) {
                                            set<int> const& s1 = eq1[j1];
                                            set<int> const& s2 = eq2[j2];
                                            auto ss(s1);
                                            ss.insert(s2.begin(), s2.end());
                                            int cs = getT(ss, sys, T, TT);
                                            if (cs < pw.second[i]) continue;
                                            auto ss1(ss);
                                            for (auto x : s1) if (ss1.count(x) != 0) ss1.erase(x);
                                            auto ss2(ss);
                                            for (auto x : s2) if (ss2.count(x) != 0) ss2.erase(x);
                                            it1->second[j1] = min(it1->second[j1], max(0, pw.second[i] - getT(ss1, sys, T, TT)));
                                            it2->second[j2] = min(it2->second[j2], max(0, pw.second[i] - getT(ss2, sys, T, TT)));
                                        }
                                    }
                                }
                            }
                        }
                    }
                    map<int, vector<int>>().swap(waiting);
                    malloc_trim(0);
                }


                // add constraints for cube variables
                for (auto x : vcube) {
                    if (setroot.count(x) == 0) {
                        GRBLinExpr e = 0;
                        for (auto y : rsys2[x]) e += X[mapVars.at(make_pair(y, x))];
                        model->addConstr(e >= 1);
                    }
                }

                // add constraints for 0s
                for (auto x : v0) {
                    for (auto y : rsys2[x]) model->addConstr(X[mapVars.at(make_pair(y, x))] == 0);
                }

                cout << endl << endl << endl << "*******************" << endl << "    solution " << ++cptt - deletedModels << "    " << endl << "*******************" << endl;
                cout << "bound: " << getT(setroot, sys, T, TT) << endl;
                cout << "setroot: " << setroot.size() << endl;
                cout << "vcube: " << vcube.size() << endl;



                // extra constraints from 1s
                model->addConstr(X[mapVars.at(make_pair(288, 285))] == 0);
                model->addConstr(X[mapVars.at(make_pair(288, 286))] == 0);
                model->addConstr(X[mapVars.at(make_pair(288, 287))] == 0);
                // root
                for (auto x : setroot) {
                    if (x >= 288) {
                        GRBLinExpr e = 0;
                        for (auto const& s : sys2.at(x)) {
                            e += X[mapVars.at(make_pair(x, *s.begin()))];
                        }
                        model->addConstr(e == 1);
                    }
                }

                int flag = 0, conCount = 0;
                for (auto const& p : sys2) {
                    GRBLinExpr e_out = 0;
                    for (auto const& s : p.second) {
                        e_out += X[mapVars[make_pair(p.first, *s.begin())]];
                    }
                    // at most one child
                    model->addConstr(e_out <= 1);


                    if (setroot.count(p.first) == 0) {
                        GRBLinExpr e_in = 0;
                        for (auto x : rsys2[p.first]) e_in += X[mapVars[make_pair(x, p.first)]];
                        model->addConstr(e_in >= e_out); // each child has at least one father
                        model->addConstr(e_in <= rsys2[p.first].size() * e_out); // each father has one child
                        if (PATTERN2) {
                            if (p.first >= 288 && sys2.find(p.first + 288) != sys2.end() && sys2.find(p.first + 2 * 288) != sys2.end()) {
                                vector<int> tmp;
                                tmp.emplace_back(*sys2[p.first][3].begin());
                                tmp.emplace_back(*sys2[p.first][3].rbegin());
                                tmp.emplace_back(*sys2[p.first + 2 * 288][3].begin());
                                tmp.emplace_back(*sys2[p.first + 2 * 288][3].rbegin());
                                if (none_of(tmp.begin(), tmp.end(), [&v0](auto x) {return v0.count(x) != 0; })) {
                                    if (all_of(tmp.begin(), tmp.end(), [](auto x) {return x < 288; })) { // add constraints at terminal nodes
                                        auto x03 = mapVars.at(make_pair(p.first, *sys2[p.first][3].begin()));
                                        auto x12 = mapVars.at(make_pair(p.first + 288, *sys2[p.first + 288][2].begin()));
                                        auto x13 = mapVars.at(make_pair(p.first + 288, *sys2[p.first + 288][3].begin()));
                                        auto x23 = mapVars.at(make_pair(p.first + 2 * 288, *sys2[p.first + 2 * 288][3].begin()));
                                        // remove pattern "3 consecutives"
                                        if (sys2.find(p.first + 3 * 288) != sys2.end()) {
                                            auto x33 = mapVars.at(make_pair(p.first + 3 * 288, *sys2[p.first + 3 * 288][3].begin()));
                                            model->addConstr(2 * X[x03] + 2 * X[x23] + X[x12] + 2 * X[x13] - X[x33] <= 4);
                                        }
                                        else
                                            model->addConstr(X[x12] + X[x13] + X[x03] + X[x23] <= 2);
                                        conCount++;
                                    }
                                }
                            }
                        }
                        if (NEWPATTERN) {
                            if (p.first >= 288 && sys2[p.first].size() < 5 && sys2.find(p.first - 288) != sys2.end() && sys2.find(p.first - 2 * 288) != sys2.end() && sys2.find(p.first - 3 * 288) != sys2.end()) {
                                vector<int> tmp;
                                tmp.emplace_back(*sys2[p.first][2].begin());
                                tmp.emplace_back(*sys2[p.first - 2 * 288][3].begin());
                                tmp.emplace_back(*sys2[p.first - 2 * 288][3].rbegin());
                                tmp.emplace_back(*sys2[p.first - 3 * 288][3].rbegin());
                                if (none_of(tmp.begin(), tmp.end(), [&v0](auto x) {return v0.count(x) != 0; })) {
                                    if (all_of(tmp.begin(), tmp.end(), [](auto x) {return x < 288; })) {
                                        auto x11 = mapVars.at(make_pair(p.first, *sys2[p.first][2].begin()));
                                        auto x22 = mapVars.at(make_pair(p.first - 1 * 288, *sys2[p.first - 1 * 288][2].begin()));
                                        auto x32 = mapVars.at(make_pair(p.first - 2 * 288, *sys2[p.first - 2 * 288][3].begin()));
                                        auto x33 = mapVars.at(make_pair(p.first - 2 * 288, *sys2[p.first - 2 * 288][2].begin()));
                                        auto x43 = mapVars.at(make_pair(p.first - 3 * 288, *sys2[p.first - 3 * 288][3].begin()));
                                        auto x44 = mapVars.at(make_pair(p.first - 3 * 288, *sys2[p.first - 3 * 288][2].begin()));
                                        model->addConstr(X[x11] + X[x32] + X[x33] + X[x43] <= 2);
                                        conCount++;
                                    }
                                }
                            }
                        }
                        if (PATTERNA) {
                            if (p.first < 288 || sys2[p.first].size() > 4) {
                                continue;
                            }
                            vector<int> roots;
                            roots.emplace_back(p.first);
                            int reg = p.first % 288;
                            switch (reg)
                            {
                            case 0:
                            {
                                roots.emplace_back(p.first - 43 * 288);
                                roots.emplace_back(p.first - 44 * 288);
                            }
                            break;
                            case 93:
                            {
                                roots.emplace_back(p.first - 25 * 288);
                                roots.emplace_back(p.first - 26 * 288);
                            }
                            break;
                            case 177:
                            {
                                roots.emplace_back(p.first - 13 * 288);
                                roots.emplace_back(p.first - 14 * 288);
                            }
                            break;
                            default:
                                break;
                            }
                            if (all_of(roots.begin(), roots.end(), [&sys2](auto x) {return sys2.find(x) != sys2.end(); })) {
                                vector<int> tmp;
                                tmp.emplace_back(*sys2[roots[0]][3].begin());
                                tmp.emplace_back(*sys2[roots[0]][3].rbegin());
                                tmp.emplace_back(*sys2[roots[2]][3].begin());
                                tmp.emplace_back(*sys2[roots[2]][3].rbegin());
                                if (none_of(tmp.begin(), tmp.end(), [&v0](auto x) {return v0.count(x) != 0; })) {
                                    if (all_of(tmp.begin(), tmp.end(), [](auto x) {return x < 288; })) {
                                        auto x11 = mapVars.at(make_pair(roots[0], *sys2[roots[0]][3].begin()));
                                        auto x21 = mapVars.at(make_pair(roots[1], *sys2[roots[1]][1].begin()));
                                        auto x22 = mapVars.at(make_pair(roots[1], *sys2[roots[1]][2].begin()));
                                        auto x32 = mapVars.at(make_pair(roots[2], *sys2[roots[2]][3].begin()));
                                        model->addConstr(X[x11] + X[x21] + X[x22] + X[x32] <= 2);
                                        conCount++;
                                    }
                                }
                            }
                        }
                        if (PATTERNB) {
                            if (p.first < 288 || sys2[p.first].size() > 4) {
                                continue;
                            }
                            vector<int> roots;
                            roots.emplace_back(p.first);
                            int reg = p.first % 288;
                            switch (reg)
                            {
                            case 0:
                            {
                                roots.emplace_back(p.first - 43 * 288);
                                roots.emplace_back(p.first - 45 * 288);
                            }
                            break;
                            case 93:
                            {
                                roots.emplace_back(p.first - 25 * 288);
                                roots.emplace_back(p.first - 27 * 288);
                            }
                            break;
                            case 177:
                            {
                                roots.emplace_back(p.first - 13 * 288);
                                roots.emplace_back(p.first - 15 * 288);
                            }
                            break;
                            default:
                                break;
                            }
                            if (all_of(roots.begin(), roots.end(), [&sys2](auto x) {return sys2.find(x) != sys2.end(); })) {
                                vector<int> tmp;
                                tmp.emplace_back(*sys2[roots[0]][3].begin());
                                tmp.emplace_back(*sys2[roots[0]][3].rbegin());
                                tmp.emplace_back(*sys2[roots[2]][3].begin());
                                tmp.emplace_back(*sys2[roots[2]][3].rbegin());
                                if (none_of(tmp.begin(), tmp.end(), [&v0](auto x) {return v0.count(x) != 0; })) {
                                    if (all_of(tmp.begin(), tmp.end(), [](auto x) {return x < 288; })) {
                                        auto x11 = mapVars.at(make_pair(roots[0], *sys2[roots[0]][3].begin()));
                                        auto x21 = mapVars.at(make_pair(roots[1], *sys2[roots[1]][1].begin()));
                                        auto x22 = mapVars.at(make_pair(roots[1], *sys2[roots[1]][2].begin()));
                                        auto x32 = mapVars.at(make_pair(roots[2], *sys2[roots[2]][3].begin()));
                                        model->addConstr(X[x11] + X[x21] + X[x22] + X[x32] <= 2);
                                        conCount++;
                                    }
                                }
                            }
                        }
                        if (PATTERNC) {
                            if (p.first < 288 || sys2[p.first].size() > 4) {
                                continue;
                            }
                            vector<int> roots;
                            roots.emplace_back(p.first);
                            int reg = p.first % 288;
                            switch (reg)
                            {
                            case 0:
                            {
                                roots.emplace_back(p.first - 44 * 288);
                                roots.emplace_back(p.first - 45 * 288);
                            }
                            break;
                            case 93:
                            {
                                roots.emplace_back(p.first - 26 * 288);
                                roots.emplace_back(p.first - 27 * 288);
                            }
                            break;
                            case 177:
                            {
                                roots.emplace_back(p.first - 14 * 288);
                                roots.emplace_back(p.first - 15 * 288);
                            }
                            break;
                            default:
                                break;
                            }
                            if (all_of(roots.begin(), roots.end(), [&sys2](auto x) {return sys2.find(x) != sys2.end(); })) {
                                vector<int> tmp;
                                tmp.emplace_back(*sys2[roots[0]][3].begin());
                                tmp.emplace_back(*sys2[roots[0]][3].rbegin());
                                tmp.emplace_back(*sys2[roots[2]][3].begin());
                                tmp.emplace_back(*sys2[roots[2]][3].rbegin());
                                if (none_of(tmp.begin(), tmp.end(), [&v0](auto x) {return v0.count(x) != 0; })) {
                                    if (all_of(tmp.begin(), tmp.end(), [](auto x) {return x < 288; })) {
                                        auto x11 = mapVars.at(make_pair(roots[0], *sys2[roots[0]][3].begin()));
                                        auto x21 = mapVars.at(make_pair(roots[1], *sys2[roots[1]][1].begin()));
                                        auto x22 = mapVars.at(make_pair(roots[1], *sys2[roots[1]][2].begin()));
                                        auto x32 = mapVars.at(make_pair(roots[2], *sys2[roots[2]][3].begin()));
                                        model->addConstr(X[x11] + X[x21] + X[x22] + X[x32] <= 2);
                                        conCount++;
                                    }
                                }
                            }
                        }
                        if (PATTERND) {
                            if (p.first < 288 || sys2[p.first].size() > 4) {
                                continue;
                            }
                            vector<int> roots;
                            roots.emplace_back(p.first);
                            int reg = p.first % 288;
                            switch (reg)
                            {
                            case 0:
                            {
                                roots.emplace_back(p.first - 44 * 288);
                                roots.emplace_back(p.first - 46 * 288);
                            }
                            break;
                            case 93:
                            {
                                roots.emplace_back(p.first - 26 * 288);
                                roots.emplace_back(p.first - 28 * 288);
                            }
                            break;
                            case 177:
                            {
                                roots.emplace_back(p.first - 14 * 288);
                                roots.emplace_back(p.first - 16 * 288);
                            }
                            break;
                            default:
                                break;
                            }
                            if (all_of(roots.begin(), roots.end(), [&sys2](auto x) {return sys2.find(x) != sys2.end(); })) {
                                vector<int> tmp;
                                tmp.emplace_back(*sys2[roots[0]][3].begin());
                                tmp.emplace_back(*sys2[roots[0]][3].rbegin());
                                tmp.emplace_back(*sys2[roots[2]][3].begin());
                                tmp.emplace_back(*sys2[roots[2]][3].rbegin());
                                if (none_of(tmp.begin(), tmp.end(), [&v0](auto x) {return v0.count(x) != 0; })) {
                                    if (all_of(tmp.begin(), tmp.end(), [](auto x) {return x < 288; })) {
                                        auto x11 = mapVars.at(make_pair(roots[0], *sys2[roots[0]][3].begin()));
                                        auto x21 = mapVars.at(make_pair(roots[1], *sys2[roots[1]][1].begin()));
                                        auto x22 = mapVars.at(make_pair(roots[1], *sys2[roots[1]][2].begin()));
                                        auto x32 = mapVars.at(make_pair(roots[2], *sys2[roots[2]][3].begin()));
                                        model->addConstr(X[x11] + X[x21] + X[x22] + X[x32] <= 2);
                                        conCount++;
                                    }
                                }
                            }
                        }
                        if (PATTERNE) {
                            if (p.first < 288 || sys2[p.first].size() > 4) {
                                continue;
                            }
                            vector<int> roots;
                            roots.emplace_back(p.first);
                            int reg = p.first % 288;
                            switch (reg)
                            {
                            case 0:
                            {
                                roots.emplace_back(p.first - 45 * 288);
                                roots.emplace_back(p.first - 46 * 288);
                            }
                            break;
                            case 93:
                            {
                                roots.emplace_back(p.first - 27 * 288);
                                roots.emplace_back(p.first - 28 * 288);
                            }
                            break;
                            case 177:
                            {
                                roots.emplace_back(p.first - 15 * 288);
                                roots.emplace_back(p.first - 16 * 288);
                            }
                            break;
                            default:
                                break;
                            }
                            if (all_of(roots.begin(), roots.end(), [&sys2](auto x) {return sys2.find(x) != sys2.end(); })) {
                                vector<int> tmp;
                                tmp.emplace_back(*sys2[roots[0]][2].begin());
                                tmp.emplace_back(*sys2[roots[2]][3].begin());
                                tmp.emplace_back(*sys2[roots[2]][3].rbegin());
                                if (none_of(tmp.begin(), tmp.end(), [&v0](auto x) {return v0.count(x) != 0; })) {
                                    if (all_of(tmp.begin(), tmp.end(), [](auto x) {return x < 288; })) {
                                        auto x11 = mapVars.at(make_pair(roots[0], *sys2[roots[0]][2].begin()));
                                        auto x21 = mapVars.at(make_pair(roots[1], *sys2[roots[1]][1].begin()));
                                        auto x22 = mapVars.at(make_pair(roots[1], *sys2[roots[1]][2].begin()));
                                        auto x32 = mapVars.at(make_pair(roots[2], *sys2[roots[2]][3].begin()));
                                        model->addConstr(X[x11] + X[x21] + X[x22] + X[x32] <= 2);
                                        conCount++;
                                    }
                                }
                            }
                        }
                        if (PATTERNF) {
                            if (p.first < 288 || sys2[p.first].size() > 4) {
                                continue;
                            }
                            vector<int> roots;
                            roots.emplace_back(p.first);
                            int reg = p.first % 288;
                            switch (reg)
                            {
                            case 0:
                            {
                                roots.emplace_back(p.first - 45 * 288);
                                roots.emplace_back(p.first - 47 * 288);
                            }
                            break;
                            case 93:
                            {
                                roots.emplace_back(p.first - 27 * 288);
                                roots.emplace_back(p.first - 29 * 288);
                            }
                            break;
                            case 177:
                            {
                                roots.emplace_back(p.first - 15 * 288);
                                roots.emplace_back(p.first - 17 * 288);
                            }
                            break;
                            default:
                                break;
                            }
                            if (all_of(roots.begin(), roots.end(), [&sys2](auto x) {return sys2.find(x) != sys2.end(); })) {
                                vector<int> tmp;
                                tmp.emplace_back(*sys2[roots[0]][2].begin());
                                tmp.emplace_back(*sys2[roots[2]][3].begin());
                                tmp.emplace_back(*sys2[roots[2]][3].rbegin());
                                if (none_of(tmp.begin(), tmp.end(), [&v0](auto x) {return v0.count(x) != 0; })) {
                                    if (all_of(tmp.begin(), tmp.end(), [](auto x) {return x < 288; })) {
                                        auto x11 = mapVars.at(make_pair(roots[0], *sys2[roots[0]][2].begin()));
                                        auto x21 = mapVars.at(make_pair(roots[1], *sys2[roots[1]][1].begin()));
                                        auto x22 = mapVars.at(make_pair(roots[1], *sys2[roots[1]][2].begin()));
                                        auto x32 = mapVars.at(make_pair(roots[2], *sys2[roots[2]][3].begin()));
                                        model->addConstr(X[x11] + X[x21] + X[x22] + X[x32] <= 2);
                                        conCount++;
                                    }
                                }
                            }
                        }
                        if (PATTERNG) {
                            if (p.first < 288 || sys2[p.first].size() > 4) {
                                continue;
                            }
                            vector<int> roots;
                            roots.emplace_back(p.first);
                            int reg = p.first % 288;
                            switch (reg)
                            {
                            case 0:
                            {
                                roots.emplace_back(p.first - 1 * 288);
                                roots.emplace_back(p.first - 43 * 288);
                            }
                            break;
                            case 93:
                            {
                                roots.emplace_back(p.first - 1 * 288);
                                roots.emplace_back(p.first - 25 * 288);
                            }
                            break;
                            case 177:
                            {
                                roots.emplace_back(p.first - 1 * 288);
                                roots.emplace_back(p.first - 13 * 288);
                            }
                            break;
                            default:
                                break;
                            }
                            if (all_of(roots.begin(), roots.end(), [&sys2](auto x) {return sys2.find(x) != sys2.end(); })) {
                                vector<int> tmp;
                                tmp.emplace_back(*sys2[roots[0]][2].begin());
                                tmp.emplace_back(*sys2[roots[0]][3].begin());
                                tmp.emplace_back(*sys2[roots[0]][3].rbegin());
                                if (none_of(tmp.begin(), tmp.end(), [&v0](auto x) {return v0.count(x) != 0; })) {
                                    if (all_of(tmp.begin(), tmp.end(), [](auto x) {return x < 288; })) {
                                        auto x11 = mapVars.at(make_pair(roots[1], *sys2[roots[1]][3].begin()));
                                        auto x21 = mapVars.at(make_pair(roots[0], *sys2[roots[0]][2].begin()));
                                        auto x22 = mapVars.at(make_pair(roots[0], *sys2[roots[0]][3].begin()));
                                        auto x32 = mapVars.at(make_pair(roots[2], *sys2[roots[2]][1].begin()));
                                        model->addConstr(X[x11] + X[x21] + X[x22] + X[x32] <= 2);
                                        conCount++;
                                    }
                                }
                            }
                        }
                        if (PATTERNH) {
                            if (p.first < 288 || sys2[p.first].size() > 4) {
                                continue;
                            }
                            vector<int> roots;
                            roots.emplace_back(p.first);
                            int reg = p.first % 288;
                            switch (reg)
                            {
                            case 0:
                            {
                                roots.emplace_back(p.first - 43 * 288);
                                roots.emplace_back(p.first - (43 + 45) * 288);
                            }
                            break;
                            case 93:
                            {
                                roots.emplace_back(p.first - 25 * 288);
                                roots.emplace_back(p.first - (25 + 27) * 288);
                            }
                            break;
                            case 177:
                            {
                                roots.emplace_back(p.first - 13 * 288);
                                roots.emplace_back(p.first - (13 + 15) * 288);
                            }
                            break;
                            default:
                                break;
                            }
                            if (all_of(roots.begin(), roots.end(), [&sys2](auto x) {return sys2.find(x) != sys2.end(); })) {
                                vector<int> tmp;
                                tmp.emplace_back(*sys2[roots[0]][3].begin());
                                tmp.emplace_back(*sys2[roots[0]][3].rbegin());
                                tmp.emplace_back(*sys2[roots[2]][1].begin());
                                if (none_of(tmp.begin(), tmp.end(), [&v0](auto x) {return v0.count(x) != 0; })) {
                                    if (all_of(tmp.begin(), tmp.end(), [](auto x) {return x < 288; })) {
                                        auto x11 = mapVars.at(make_pair(roots[0], *sys2[roots[0]][3].begin()));
                                        auto x21 = mapVars.at(make_pair(roots[1], *sys2[roots[1]][1].begin()));
                                        auto x22 = mapVars.at(make_pair(roots[1], *sys2[roots[1]][2].begin()));
                                        auto x32 = mapVars.at(make_pair(roots[2], *sys2[roots[2]][1].begin()));
                                        model->addConstr(X[x11] + X[x21] + X[x22] + X[x32] <= 2);
                                        conCount++;
                                    }
                                }
                            }
                        }
                        if (PATTERNI) {
                            if (p.first < 288 || sys2[p.first].size() > 4) {
                                continue;
                            }
                            vector<int> roots;
                            roots.emplace_back(p.first);
                            int reg = p.first % 288;
                            switch (reg)
                            {
                            case 0:
                            {
                                roots.emplace_back(p.first - 44 * 288);
                                roots.emplace_back(p.first - (44 + 45) * 288);
                            }
                            break;
                            case 93:
                            {
                                roots.emplace_back(p.first - 26 * 288);
                                roots.emplace_back(p.first - (26 + 27) * 288);
                            }
                            break;
                            case 177:
                            {
                                roots.emplace_back(p.first - 14 * 288);
                                roots.emplace_back(p.first - (14 + 15) * 288);
                            }
                            break;
                            default:
                                break;
                            }
                            if (all_of(roots.begin(), roots.end(), [&sys2](auto x) {return sys2.find(x) != sys2.end(); })) {
                                vector<int> tmp;
                                tmp.emplace_back(*sys2[roots[0]][3].begin());
                                tmp.emplace_back(*sys2[roots[0]][3].rbegin());
                                tmp.emplace_back(*sys2[roots[2]][1].begin());
                                if (none_of(tmp.begin(), tmp.end(), [&v0](auto x) {return v0.count(x) != 0; })) {
                                    if (all_of(tmp.begin(), tmp.end(), [](auto x) {return x < 288; })) {
                                        auto x11 = mapVars.at(make_pair(roots[0], *sys2[roots[0]][3].begin()));
                                        auto x21 = mapVars.at(make_pair(roots[1], *sys2[roots[1]][1].begin()));
                                        auto x22 = mapVars.at(make_pair(roots[1], *sys2[roots[1]][2].begin()));
                                        auto x32 = mapVars.at(make_pair(roots[2], *sys2[roots[2]][1].begin()));
                                        model->addConstr(X[x11] + X[x21] + X[x22] + X[x32] <= 2);
                                        conCount++;
                                    }
                                }
                            }
                        }
                        if (PATTERNJ) {
                            if (p.first < 288 || sys2[p.first].size() > 4) {
                                continue;
                            }
                            vector<int> roots;
                            roots.emplace_back(p.first);
                            int reg = p.first % 288;
                            switch (reg)
                            {
                            case 0:
                            {
                                roots.emplace_back(p.first - 45 * 288);
                                roots.emplace_back(p.first - 2 * 45 * 288);
                            }
                            break;
                            case 93:
                            {
                                roots.emplace_back(p.first - 27 * 288);
                                roots.emplace_back(p.first - 2 * 27 * 288);
                            }
                            break;
                            case 177:
                            {
                                roots.emplace_back(p.first - 15 * 288);
                                roots.emplace_back(p.first - 2 * 15 * 288);
                            }
                            break;
                            default:
                                break;
                            }
                            if (all_of(roots.begin(), roots.end(), [&sys2](auto x) {return sys2.find(x) != sys2.end(); })) {
                                vector<int> tmp;
                                tmp.emplace_back(*sys2[roots[0]][2].begin());
                                tmp.emplace_back(*sys2[roots[1]][2].begin());
                                if (none_of(tmp.begin(), tmp.end(), [&v0](auto x) {return v0.count(x) != 0; })) {
                                    if (all_of(tmp.begin(), tmp.end(), [](auto x) {return x < 288; })) {
                                        auto x11 = mapVars.at(make_pair(roots[0], *sys2[roots[0]][2].begin()));
                                        auto x21 = mapVars.at(make_pair(roots[1], *sys2[roots[1]][1].begin()));
                                        auto x22 = mapVars.at(make_pair(roots[1], *sys2[roots[1]][2].begin()));
                                        auto x32 = mapVars.at(make_pair(roots[2], *sys2[roots[2]][1].begin()));
                                        model->addConstr(X[x11] + X[x21] + X[x22] + X[x32] <= 2);
                                        conCount++;
                                    }
                                }
                            }
                        }
                    }
                }


                // model.read("tune.prm");
                model->set(GRB_IntParam_RINS, 0);
                model->set(GRB_IntParam_MIPFocus, 3);
                model->set(GRB_IntParam_VarBranch, 2);
                // model.set(GRB_IntParam_OutputFlag, 0);
                model->set(GRB_DoubleParam_MIPGap, GRB_INFINITY);
                model->set(GRB_IntParam_PoolSearchMode, 2);
                model->set(GRB_IntParam_PoolSolutions, 2000000000);
                // model.set("TimeLimit", timelimit);
                {
                    if (vsetroot.size() - deletedModels < 5000) {
                        model->set("TimeLimit", "150");
                    }
                    else if (vsetroot.size() - deletedModels < 10000) {
                        model->set("TimeLimit", "300");
                    }
                    else if (vsetroot.size() - deletedModels < 20000) {
                        model->set("TimeLimit", "600");
                    }
                    else if (vsetroot.size() - deletedModels < 30000) {
                        model->set("TimeLimit", "1200");
                    }
                    else if (vsetroot.size() - deletedModels < 40000) {
                        model->set("TimeLimit", "2400");
                    }
                    else{
                        model->set("TimeLimit", "3600");
                    }
                }
                /*if (SAVE_MILP) {
                    string filename = "round840_root";
                    char rootStr[10];
                    sprintf(rootStr, "%d", *setroot.begin());
                    model->write((filename.append(rootStr)).append(".lp"));
                }*/
                printf("Solving...\n");
                model->optimize();
                if (model->get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
                    printf("This solution is time out!\n");
                    printf("Deviding...\n");
                    if (setroot.size() == 1) {
                        for (auto const& s : sys.at(*setroot.begin())) {
                            if (*s.begin() == *s.rbegin())
                            {
                                if (getT(*s.begin(), sys, T, TT) >= vcube.size()) {
                                    vsetroot.emplace_back(s);
                                }
                            }
                            if (*s.begin() != *s.rbegin())
                            {
                                if (getT(*s.begin(), sys, T, TT) + getT(*s.rbegin(), sys, T, TT) >= vcube.size()) {
                                    vsetroot.emplace_back(s);
                                }
                            }
                        }
                    }
                    else if (setroot.size() <= 2) {
                        int rootSize = setroot.size();
                        int* setrootarray = (int*)malloc(rootSize * sizeof(int));
                        set<int>** subnodearray = (set<int>**)malloc(4 * sizeof(set<int>*));
                        int index = 0;
                        while (!setroot.empty()) {
                            setrootarray[index] = *setroot.begin();
                            setroot.erase(setrootarray[index++]);
                        }
                        for (int i = 0; i < pow(4, rootSize); i++) {
                            set<int> subRoot;
                            int highnum = 0;
                            for (int j = 0; j < rootSize; j++) {
                                int k = 0;
                                for (auto& subnode : sys2.at(setrootarray[j])) {
                                    subnodearray[k++] = &subnode;
                                }
                                int dim = (int)((i - highnum) / pow(4, rootSize - j - 1));
                                if (dim < 3) {
                                    subRoot.insert(*subnodearray[dim]->begin());
                                }
                                else {
                                    subRoot.insert(*subnodearray[dim]->begin());
                                    subRoot.insert(*subnodearray[dim]->rbegin());
                                }
                                highnum += (int)(dim * pow(4, rootSize - j - 1));
                            }
                            if (getT(subRoot, sys, T, TT) >= vcube.size()) {
                                vsetroot.emplace_back(subRoot);
                            }
                            set<int>().swap(subRoot);
                            malloc_trim(0);
                        }
                    }
                    else {
                        int rootSize = setroot.size();
                        int* setrootarray = (int*)malloc(rootSize * sizeof(int));
                        // set<int>** subnodearray = (set<int>**)malloc(4 * sizeof(set<int>*));
                        int index = 0;
                        while (!setroot.empty()) {
                            setrootarray[index] = *setroot.begin();
                            setroot.erase(setrootarray[index++]);
                        }
                        pthread_mutex_t mutex;
                        pthread_barrier_t barrier;
                        pthread_mutex_init(&mutex, NULL);
                        pthread_barrier_init(&barrier, NULL, MULTITHREAD);
                        pthread_t* threads = (pthread_t*)malloc(sizeof(pthread_t) * MULTITHREAD);
                        for (int i = 0; i < MULTITHREAD; i++) {
                            args4MultiThread* args = (args4MultiThread*)malloc(sizeof(args4MultiThread));
                            args->rootSize = rootSize;
                            args->setrootarray = setrootarray;
                            // args->subnodearray = subnodearray;
                            args->UID = i;
                            args->vcubesize = vcube.size();
                            args->vsetroot = &vsetroot;
                            args->T = &T;
                            args->TT = &TT;
                            args->sys = &sys;
                            args->sys2 = &sys2;
                            args->mutex = &mutex;
                            args->barrier = &barrier;
                            pthread_create(threads + i, NULL, devideMILP, (void*)(args));
                        }
                        for (int i = 0; i < MULTITHREAD; i++) {
                            pthread_join(threads[i], NULL);
                        }
                        pthread_mutex_destroy(&mutex);
                        pthread_barrier_destroy(&barrier);
                    }
                    deletedModels++;
                    delete model;
                    map<int, vector<set<int>>>().swap(sys2);
                    map<int, vector<int>>().swap(rsys2);
                    map<pair<int, int>, int>().swap(mapVars);
                    malloc_trim(0);
                    continue;
                }

                nSolutions += model->get(GRB_IntAttr_SolCount);
                cout << "\r" << "cptt: " << cptt - deletedModels << "/" << vsetroot.size() - deletedModels << " -- " << "solutions: " << nSolutions << endl;

                //save superpoly
                int solCount = model->get(GRB_IntAttr_SolCount);
                if (solCount == 0) {
                    delete model;
                    map<int, vector<set<int>>>().swap(sys2);
                    map<int, vector<int>>().swap(rsys2);
                    map<pair<int, int>, int>().swap(mapVars);
                    malloc_trim(0);
                    continue;
                }
                if (solCount >= 2000000000) {
                    cerr << "Number of solutions is too large" << endl;
                    exit(0);
                }
                if (PRINT_GRAPH) {
                    for (int i = 0; i < 10; i++) {
                        int numvars = model->get(GRB_IntAttr_NumVars);
                        model->set(GRB_IntParam_SolutionNumber, i);
                        set<pair<int, int>> edges;
                        for (int j = 0; j < numvars; j++) {
                            auto v = &X[j];
                            if (v->get(GRB_DoubleAttr_Xn) > 0.5) {
                                int k = 0;
                                for (; k < X.size(); k++) {
                                    if (&(X[k]) == v)
                                    {
                                        break;
                                    }
                                }
                                for (auto x : mapVars) {
                                    if (x.second == k)
                                    {
                                        edges.insert(x.first);
                                    }
                                }
                            }
                        }
                        for (auto edge : edges) {
                            cout << "<" << edge.first << ", " << edge.second << ">" << endl;
                        }
                        cout << endl;
                    }
                }
                printf("Couting...\n");
                for (int i = 0; i < solCount; i++) {
                    model->set(GRB_IntParam_SolutionNumber, i);
                    bitset<288> tmp;
                    for (int x = 0; x < 288; x++) {
                        tmp[x] = 0;
                        for (auto y : rsys.at(x)) {
                            if (mapVars.find(make_pair(y, x)) != mapVars.end()) {
                                if (X[mapVars.at(make_pair(y, x))].get(GRB_DoubleAttr_Xn) > 0.5) {
                                    tmp[x] = 1;
                                }
                            }
                        }
                    }
                    countingBox[tmp]++;
                }
                delete model;
                map<int, vector<set<int>>>().swap(sys2);
                map<int, vector<int>>().swap(rsys2);
                map<pair<int, int>, int>().swap(mapVars);
                malloc_trim(0);
            }
        }
        cout << "SUPERPOLY" << endl;
        display_result_anf(vcube, countingBox);
        cout << endl;
    }
    catch (GRBException e) {
        cout << "Error  number: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
    catch (...) {
        cout << "Error  during  optimization" << endl;
    }
}

int main(int argc, char const* argv[]) {
    int R = stoi(argv[1]);
    set<int> vcube, v0;
    int backRound = 0;
    switch (R) {
    case 842:
    {
        for (int i = 0; i < 80; ++i) {
            if (i == 18 || i == 34) v0.emplace(i);
            else vcube.emplace(i);
        }
        backRound = 320;
    }
    break;
    case 843:
    {
        vcube = set<int>({  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                           21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 34, 36, 38, 40, 42, 45, 47, 49, 51,
                           53, 55, 57, 60, 62, 64, 66, 68, 70, 72, 77, 75, 79
            });
        for (int i = 0; i < 80; ++i) {
            if (vcube.count(i) == 0) v0.emplace(i);
        }
        backRound = 320;
    }
    break;
    case 845:
    {
        vcube = set<int>({  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, // I3
                           21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 34, 36, 38, 40, 42, 45, 47, 49, 51, 53, 
                           55, 57, 60, 62, 64, 66, 68, 70, 72, 77, 75, 79 });
        for (int i = 0; i < 80; ++i) {
            if (vcube.count(i) == 0) v0.emplace(i);
        }
        backRound = 320;
    }
    break;
    default:
    {
    }
    }
    auto start = high_resolution_clock::now();
    milpTrivium(R, vcube, v0, backRound);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    cout << "Time: " << duration.count() / 1000000 << "s" << endl;

    return 0;
}
