#include<bits/stdc++.h>
#include<unordered_map>
#include "cdfppf.h"
using namespace std;
typedef pair<double, double> DD;
typedef pair<double, int> DI;
typedef pair<int, vector<DD>> IV;
typedef pair<int, vector<int>> IVI;
typedef pair<int, int> II;
typedef pair<int, double> ID;
const int MAX_V = 435676, MAX_M = 1057266;
const int populationsize = 10, generations = 20;
int N, M;
double optw = DBL_MAX;
vector<DD> coord;
double coord2trav = DBL_MAX;

struct edge{
    int to, id;
    double mean, var;
    double realizew[10];
    edge(int _to, int _id, double _mean, double _var){
        to = _to;
        id = _id;
        mean = _mean;
        var = _var;
    }
};

struct path{
    int lastid, to;
    double mean, var;
    path(int _id, int _to, double _mean, double _var){
        lastid = _id;
        to = _to;
        mean = _mean;
        var = _var;
    }
};

vector<edge> adj[MAX_V];
unordered_map<int,double> cov[MAX_M];

double fetchFinvValue(double alpha){
    alpha = (alpha < 0.001) ? 0.001 : alpha;
    return ppf[(int)(alpha * 1000) - 1];
}

double fetchTrav(double mean, double var, double alpha){
    return mean + fetchFinvValue(alpha) * sqrt(var);
}

vector<vector<int>> popu;

double fitness(int s, vector<int> &path, double alpha){
    int v = s;
    double tub[10] = {0};
    for(int i = 0; i < path.size(); i++){
        edge e = adj[v][path[i]];
        for (int j = 0; j < 10;j++)
            tub[j] += e.realizew[j];
        v = e.to;
    }
    sort(tub, tub + 10);
    return tub[int(alpha * 10)];
}

double travtime(int s, vector<int> &path, double alpha){
    int v = s;
    double mean = 0, var = 0;
    for(int i = 0; i < path.size(); i++){
        edge e = adj[v][path[i]];
        v = e.to;
        mean += e.mean;
        var += e.var;
        if(i!=path.size()-1)
            var += 2 * cov[e.id][adj[v][path[i + 1]].id];
    }
    return fetchTrav(mean, var, alpha);
}

typedef pair<double, int> P;
int pre[MAX_V], pre_i[MAX_V];
double dist[MAX_V];
vector<int> shortestpath(int s, int t, int jwei){
    priority_queue<P, vector<P>, greater<P> > que;
    for (int i = 0; i <= N;i++)
        dist[i] = DBL_MAX;
    memset(pre, -1, sizeof(pre));
    memset(pre_i, -1, sizeof(pre_i));
    dist[s] = 0;
    que.push(P(0, s));
    vector<int> path;
    while (!que.empty()){
        P p = que.top();
        que.pop();
        int v = p.second;
        if (dist[v] < p.first)
            continue;
        if (v == t){
            for (; v != s; v = pre[v])
                path.push_back(pre_i[v]);
            reverse(path.begin(), path.end());
            return path;
        }
        for (int i = 0; i < adj[v].size(); i++) {
            edge e = adj[v][i];
            if (dist[e.to] > dist[v] + e.realizew[jwei]){
                dist[e.to] = dist[v] + e.realizew[jwei];
                pre[e.to] = v;
                pre_i[e.to] = i;
                que.push(P(dist[e.to], e.to));
            }
        }
    }
    return path;
}

int vflag[MAX_V];
void crossover(int s, vector<int> &P1,vector<int> &P2){
    memset(vflag, 0, sizeof(vflag));
    int v = s;
    for(int i = 0; i < P1.size(); i++){
        edge e = adj[v][P1[i]];
        v = e.to;
        vflag[v] = i;
        //printf("%d ", v);
    }
    vflag[v] = 0;
    v = s;
    vector<II> commonnode;
    for(int i = 0; i < P2.size(); i++){
        edge e = adj[v][P2[i]];
        v = e.to;
        //printf("%d ", v);
        if(vflag[v])
            commonnode.push_back(II(v, i));
    }
    if (commonnode.size() == 0)
        return;
    int randind = rand() % commonnode.size();
    int P2ind = commonnode[randind].second;
    int P1ind = vflag[commonnode[randind].first];
    //printf("%d %d %d\n", P1ind, P2ind, commonnode[randind].first);
    vector<int> P1off, P2off;
    for (int i = 0; i < P1ind + 1; i++){
        edge e = adj[v][P1[i]];
        v = e.to;
        P1off.push_back(P1[i]);
    }
    for (int i = P2ind + 1; i < P2.size(); i++){
        P1off.push_back(P2[i]);
    }
    for (int i = 0; i < P2ind + 1; i++){
        edge e = adj[v][P2[i]];
        v = e.to;
        P2off.push_back(P2[i]);
    }
    for (int i = P1ind + 1; i < P1.size(); i++){
        P2off.push_back(P1[i]);
    }
    popu.push_back(P1off);
    popu.push_back(P2off);
}

void offsprings(int s){
    while(popu.size()!=populationsize)
        popu.pop_back();
    int len2popu = popu.size() / 2;
    for (int i = 0; i < len2popu; i++)
        popu.pop_back();
    int lenpopu = popu.size();
    for (int i = 0; i < lenpopu;i++){
        for (int j = i + 1; j < lenpopu; j++){
            //printf("%d %d\n", i, j);
            crossover(s, popu[i], popu[j]);
        }
    }
}

double SMOGAQuery(int s, int t, double alpha){
    s--;
    t--;
    popu.clear();
    if (s == t)
        return 0;
    for (int i = 0; i < populationsize;i++)
        popu.push_back(shortestpath(s, t, i));
    for (int i = 0; i < generations;i++){
        vector<P> sortedP;
        for (int j = 0; j < popu.size();j++){
            sortedP.push_back(P(fitness(s, popu[j], alpha), j));
        }
        sort(sortedP.begin(), sortedP.end());
        offsprings(s);
    }
    double minp = DBL_MAX;
    int ind = -1;
    for (int j = 0; j < popu.size();j++){
        if(fitness(s, popu[j], alpha)<minp){
            ind = j;
        }
    }
    return travtime(s, popu[ind], alpha);
}

int main(int argc , char * argv[]){
    string sfile, sq, st;
    FILE *fp_query, *fp_networkvar, *fp_networkmu, *fp_networkcov, *fp_networkcoord;
    if (argc > 1)
        sfile = string(argv[1]);
    else
        sfile = string("NY");
    if (argc > 2)
        sq = string(argv[2]);
    else
        sq = string("");
    if (argc > 3){
        st = string(argv[3]);
    }
    else
        st = string("");
    double war = cdf[0];
    string prefix = string("../data/") + sfile + string("/");
    string smu = prefix + string("USA-road-t.") + sfile + (".gr");
    string svar = prefix + sq+string("USA-road-var.") + sfile + (".gr");
    string scov = prefix + st+string("USA-road-cov.") + sfile + (".gr");
    string scoord = prefix + string("USA-road-d.") + sfile + (".co");

    fp_networkmu = fopen(smu.c_str(), "r");
    fp_networkvar = fopen(svar.c_str(), "r");
    fp_networkcov = fopen(scov.c_str(), "r");
    fp_networkcoord = fopen(scoord.c_str(), "r");
    char ch, buffer[100];
    int u, v;
    double w, c;
    //travel time
    for (int i = 0; i < 4; i++){
        fgets(buffer, 90, fp_networkvar);
        fgets(buffer, 90, fp_networkmu);
        fgets(buffer, 90, fp_networkcoord);
    }
    for (int i = 0; i < 4; i++){
        fgetc(fp_networkvar);
        fgetc(fp_networkmu);
    }
    fscanf(fp_networkmu, "%d%d", &N, &M);
    fscanf(fp_networkvar, "%d%d", &N, &M);
    for (int i = 0; i < 3; i++){
        fgets(buffer, 90, fp_networkvar);
        fgets(buffer, 90, fp_networkmu);
        fgets(buffer, 90, fp_networkcoord);
    }
    for (int i = 0; i < N; i++){
        fscanf(fp_networkcoord, "%c%d%lf%lf", &ch, &u, &w, &c);
        fgets(buffer, 90, fp_networkcoord);
        coord.push_back(DD(w/1000000, c/1000000));
    }
    for (int i = 0; i < M; i++) {
        fscanf(fp_networkvar, "%c%d%d%lf", &ch, &u, &v, &w);
        fgets(buffer, 90, fp_networkvar);
        fscanf(fp_networkmu, "%c%d%d%lf", &ch, &u, &v, &c);
        fgets(buffer, 90, fp_networkmu);
        u--;
        v--;
        if(i%2==0){
            adj[u].push_back(edge(v, i / 2, c, w));
            adj[v].push_back(edge(u, i / 2, c, w));
            std::default_random_engine gen;
            std::normal_distribution<double> dis(c, sqrt(w));
            for (int j = 0; j < 10;j++){
                double tmp = max(0.0, dis(gen));
                adj[u][adj[u].size() - 1].realizew[j] = tmp;
                adj[v][adj[v].size() - 1].realizew[j] = tmp;
            }
        }
    }    


    while (fscanf(fp_networkcov, "%d%d%lf", &u, &v, &c) == 3){
        cov[u][v] = cov[v][u] = c;
    }

    std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT;

    freopen((prefix + sfile+sq+st+string("SMOGAResults")).c_str(), "w", stdout);
    for (int qi = 0; qi < 5;qi++){
        string s3 = string("../data/") + sfile + string("/") + string("dis") + to_string(qi+1);
        fp_query = fopen(s3.c_str(), "r");
        vector<pair<II, double>> queryset;
        int qs, qt;
        double qC;
        while (~fscanf(fp_query, "%d%d%lf", &qs, &qt, &qC)){
            queryset.push_back(make_pair(II(qs, qt), qC));
        }
        t1=std::chrono::high_resolution_clock::now();
        for(int i=0;i<1;i++)
        SMOGAQuery(queryset[i].first.first, queryset[i].first.second, queryset[i].second);
        t2=std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        runT= time_span.count();
        cout<<"Query Time "<<runT<<endl;
    }
    printf("\n\n");
    for (int qi = 0; qi < 5;qi++){
        string s3 = string("../data/") + sfile + string("/") + string("alpha") + to_string(qi+1);
        fp_query = fopen(s3.c_str(), "r");
        vector<pair<II, double>> queryset;
        int qs, qt;
        double qC;
        while (~fscanf(fp_query, "%d%d%lf", &qs, &qt, &qC)){
            queryset.push_back(make_pair(II(qs, qt), qC));
        }
        t1=std::chrono::high_resolution_clock::now();
        for(int i=0;i<1;i++)
        SMOGAQuery(queryset[i].first.first, queryset[i].first.second, queryset[i].second);
        t2=std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        runT= time_span.count();
        cout<<"Query Time "<<runT<<endl;
    }
}