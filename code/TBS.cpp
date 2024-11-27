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
    double elb;
    edge(int _to, int _id, double _mean, double _var){
        to = _to;
        id = _id;
        mean = _mean;
        var = _var;
        elb = max(0.0, mean - 2 * sqrt(var));
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

double prio(int v,int t){
    return 0;
}

double df[MAX_V], db[MAX_V];

typedef struct pqn{
	double Aw;
    double mean, var;
    int v, lastid;
    pqn(double _Aw, int _v,double _m, double _var, int _lastid){
        Aw = _Aw, mean = _m, v = _v, var = _var, lastid = _lastid;
    }
};
struct cmppqn{
    bool operator()(const pqn &a, const pqn &b){
        return a.Aw > b.Aw;
    }
};
double travub;
void rev(int s, int t, double alpha){
    priority_queue<pqn, vector<pqn>, cmppqn> Qf;
    for (int i = 0; i <= N;i++)
        db[i] = DBL_MAX;
    db[t] = 0;
    Qf.push(pqn(0, t, 0, 0, 0));
    while (!Qf.empty()){
        pqn p = Qf.top();
        Qf.pop();
        int v = p.v;
        if (db[v] < p.Aw)
            continue;
        if (v == s){
            travub = fetchTrav(p.mean, p.var, alpha);
            return;
        }
        for (int i = 0; i < adj[v].size(); i++) {
            edge e = adj[v][i];
            if (db[e.to] > db[v] + e.elb){
                db[e.to] = db[v] + e.elb;
                Qf.push(pqn(db[e.to], e.to, p.mean + e.mean, p.var + e.var, 0));
            }
        }
    }
    travub = DBL_MAX;
}

double TBSQuery(int s, int t, double alpha){
    if (s == t)
        return optw = 0;
    priority_queue<pqn, vector<pqn>, cmppqn> Qf;
    for (int i = 0; i <= N;i++)
        df[i] = DBL_MAX;
    df[s] = 0;
    Qf.push(pqn(0, s, 0, 0, -1));
    while (!Qf.empty()){
        pqn p = Qf.top();
        Qf.pop();
        int v = p.v;
        if (df[v] < p.Aw)
            continue;
        if (v == t){
            return optw = fetchTrav(p.mean, p.var, alpha);
        }
        for (int i = 0; i < adj[v].size(); i++) {
            edge e = adj[v][i];
            double travt, covar = 0;
            if (p.lastid != -1)
                covar = 2 * cov[p.lastid][e.id];
            travt = fetchTrav(p.mean + e.mean, p.var + e.var + covar, alpha);
            if(travt+db[e.to]>travub)
                continue;
            if (df[e.to] > travt){
                df[e.to] = travt;
                Qf.push(pqn(df[e.to], e.to, p.mean + e.mean, p.var + e.var + covar, e.id));
            }
        }
    }
    return 0;
}


int main(int argc , char * argv[]){
    string sfile, sq,st;
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
        }
    }
    
    while (fscanf(fp_networkcov, "%d%d%lf", &u, &v, &c) == 3){
        cov[u][v] = cov[v][u] = c;
    }

    std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT = 0;
	
	if (argc > 4){
        string s3 = prefix + string(argv[4]);
        fp_query = fopen(s3.c_str(), "r");
        vector<pair<II, double>> queryset;
        int qs, qt;
        double qC;
        while (~fscanf(fp_query, "%d%d%lf", &qs, &qt, &qC)){
            queryset.push_back(make_pair(II(qs-1, qt-1), qC));
        }

        vector<double> ans;
        for (int i = 0; i < queryset.size(); i++){
            rev(queryset[i].first.first, queryset[i].first.second, queryset[i].second);
            t1=std::chrono::high_resolution_clock::now();
            TBSQuery(queryset[i].first.first, queryset[i].first.second, queryset[i].second);
            t2=std::chrono::high_resolution_clock::now();
            time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
            runT += time_span.count();
            ans.push_back(optw);
        }
		
		FILE *fp_ans1 = fopen((s3 + string("RSPans_TBS")).c_str(), "w");
        for (int i = 0; i < ans.size();i++)
            fprintf(fp_ans1, "%f\n", ans[i]);
        fclose(fp_ans1);

        cout<<"Query Time "<<runT<<endl;
        return 0;
    }



    //freopen((prefix + sfile+sq+st+string("TBSResults")).c_str(), "w", stdout);
    printf("%s\n", (sfile + sq + st).c_str());
    for (int qi = 0; qi < 5;qi++){
        string s3 = string("../data/") + sfile + string("/") + string("dis") + to_string(qi+1);
        fp_query = fopen(s3.c_str(), "r");
        vector<pair<II, double>> queryset;
        int qs, qt;
        double qC;
        while (~fscanf(fp_query, "%d%d%lf", &qs, &qt, &qC)){
            queryset.push_back(make_pair(II(qs-1, qt-1), qC));
        }
        runT = 0;
        for (int i = 0; i < queryset.size(); i++){
            rev(queryset[i].first.first, queryset[i].first.second, queryset[i].second);
            t1=std::chrono::high_resolution_clock::now();
            TBSQuery(queryset[i].first.first, queryset[i].first.second, queryset[i].second);
            t2=std::chrono::high_resolution_clock::now();
            time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
            runT += time_span.count();
        }
        cout<<"Query Time "<<runT<<endl;
    }
    printf("\n");
    for (int qi = 0; qi < 5;qi++){
        string s3 = string("../data/") + sfile + string("/") + string("alpha") + to_string(qi+1);
        fp_query = fopen(s3.c_str(), "r");
        vector<pair<II, double>> queryset;
        int qs, qt;
        double qC;
        while (~fscanf(fp_query, "%d%d%lf", &qs, &qt, &qC)){
            queryset.push_back(make_pair(II(qs-1, qt-1), qC));
        }
        runT = 0;
        for (int i = 0; i < queryset.size(); i++){
            rev(queryset[i].first.first, queryset[i].first.second, queryset[i].second);
            t1=std::chrono::high_resolution_clock::now();
            TBSQuery(queryset[i].first.first, queryset[i].first.second, queryset[i].second);
            t2=std::chrono::high_resolution_clock::now();
            time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
            runT += time_span.count();
        }
        cout<<"Query Time "<<runT<<endl;
    }
    printf("\n");
}