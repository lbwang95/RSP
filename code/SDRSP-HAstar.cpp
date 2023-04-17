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
int N, M;
double optw = DBL_MAX;
vector<DD> coord;
double coord2trav = DBL_MAX;

struct edge{
    int to, id;
    double mean, var;
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

vector<int> Pe[MAX_M];
vector<int> validflagP;
vector<path> AllP;
vector<edge> adj[MAX_V];
unordered_map<int,double> cov[MAX_M];
vector<int> order;

double fetchFinvValue(double alpha){
    alpha = (alpha < 0.001) ? 0.001 : alpha;
    return ppf[(int)(alpha * 1000) - 1];
}

double fetchTrav(double mean, double var, double alpha){
    return mean + fetchFinvValue(alpha) * sqrt(var);
}

double distance(int u, int v){
    return sqrt(pow(coord[u].first - coord[v].first, 2) + pow(coord[u].second - coord[v].second, 2));
}

double h(int u,int v){
    return distance(u, v) * coord2trav;
}

struct alledge{
    int from, to;
    double mean, variance;
    alledge(int a,int b,double c,double d){
        from = a, to = b, mean = c, variance = d;
    }
};
vector<alledge> alledges;

typedef pair<double, int> P;
priority_queue<P, vector<P>, greater<P> > SE;

bool checkDom(int lasteid, int to, double smean, double svar, double alpha){
    double beta;
    if(alpha>0.5)
        beta = 0.999;
    if(alpha==0.5)
        beta = 0.5;
    if(alpha<0.5)
        beta = 0.001;
    int pos = 0;
    for (int i = 0; i < Pe[lasteid].size();i++){
        int iterpid = Pe[lasteid][i];
        double phip = fetchTrav(smean, svar, beta);
        double phiIterp = fetchTrav(AllP[iterpid].mean, AllP[iterpid].var, beta);
        if(smean>=AllP[iterpid].mean){
            if(phip>=phiIterp){
                return false;
            }
            pos += 1;
        }
    }
    vector<int> append;
    int lenPe = Pe[lasteid].size();
    for (int i = pos; i < lenPe; i++){
        int iterpid = Pe[lasteid][i];
        double phip = fetchTrav(smean, svar, beta);
        double phiIterp = fetchTrav(AllP[iterpid].mean, AllP[iterpid].var, beta);
        if(phip>phiIterp){
            append.push_back(iterpid);
        }
        else
            validflagP[iterpid] = 1;
    }
    for (int i = 0; i < lenPe - pos;i++)
        Pe[lasteid].pop_back();
    validflagP.push_back(0);
    AllP.push_back(path(lasteid, to, smean, svar));
    Pe[lasteid].push_back(AllP.size() - 1);
    for (int i = 0; i < append.size();i++)
        Pe[lasteid].push_back(append[i]);
    return true;
}

double ERSPQuery(int s, int t, double alpha){
    s--;
    t--;
    if (s == t)
        return 0;
    AllP.clear();
    while(SE.size())
        SE.pop();
    for (int i = 0; i < M;i++)
        Pe[i].clear();
    validflagP.clear();
    for (int i = 0; i < adj[s].size();i++){
        int v = adj[s][i].to, eid = adj[s][i].id;
        double trav = adj[s][i].mean, travvar = adj[s][i].var;
        //printf("(%d %d %f %f)\n", v, eid, trav, travvar);
        AllP.push_back(path(eid, v, trav, travvar));
        Pe[eid].push_back(AllP.size() - 1);
        validflagP.push_back(0);
        SE.push(P(h(v, t) + fetchTrav(trav, travvar, alpha), AllP.size() - 1));
    }
    while(!SE.empty()){
        P p = SE.top();
        SE.pop();
        int pid = p.second;
        if(validflagP[pid])
            continue;
        int v = AllP[pid].to, lasteid = AllP[pid].lastid;
        if (v == t)
            return p.first;
        for (int i = 0; i < adj[v].size();i++){
            if(adj[v][i].id==lasteid)
                continue;
            int w = adj[v][i].to;
            int thiseid = adj[v][i].id;
            double variance = AllP[pid].var + adj[v][i].var + 2 * cov[lasteid][thiseid];
            if(checkDom(thiseid, w, AllP[pid].mean + adj[v][i].mean, variance, alpha)){
                SE.push(P(h(w, t) + fetchTrav(AllP[pid].mean + adj[v][i].mean, variance, alpha), AllP.size() - 1));
            }
        }
    }

    return 0;
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
    //test or not
    if(0){
        sq = string("test");
        smu = string("../data/") + sfile + string("/t.") + sfile + (".gr");
    }
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
            alledges.push_back(alledge(u, v, c, w));
            coord2trav = min(coord2trav, c / distance(u, v));
        }
    }
    
    while (fscanf(fp_networkcov, "%d%d%lf", &u, &v, &c) == 3){
        cov[u][v] = cov[v][u] = c;
    }

    std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT;
    
	if (argc > 4){
        string s3 = prefix + string(argv[4]);
        fp_query = fopen(s3.c_str(), "r");
        vector<pair<II, double>> queryset;
        int qs, qt;
        double qC;
        while (~fscanf(fp_query, "%d%d%lf", &qs, &qt, &qC)){
            queryset.push_back(make_pair(II(qs, qt), qC));
        }

        vector<double> ans;
        t1=std::chrono::high_resolution_clock::now();
        for (int i = 0; i < queryset.size(); i++){
            ERSPQuery(queryset[i].first.first, queryset[i].first.second, queryset[i].second);
            ans.push_back(optw);
        }
        t2=std::chrono::high_resolution_clock::now();

        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        runT= time_span.count();
        cout<<"Query Time "<<runT<<endl;
        return 0;
    }


    freopen((prefix + sfile+sq+st+string("SDRSPResults")).c_str(), "w", stdout);
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
        for(int i=0;i<100;i++)
        ERSPQuery(queryset[i].first.first, queryset[i].first.second, queryset[i].second);
        t2=std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        runT= time_span.count();
        cout<<"Query Time "<<runT/100<<endl;
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
        for(int i=0;i<100;i++)
        ERSPQuery(queryset[i].first.first, queryset[i].first.second, queryset[i].second);
        t2=std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        runT= time_span.count();
        cout<<"Query Time "<<runT/100<<endl;
    }
}