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
const int MAX_V = 435676;//1070576
int N, M;
long long npathConcat, hopsize;
double optw = DBL_MAX;
int treeheight = 0, treewidth = 0, treeavgheight = 0;

double FinvValue(double alpha){//the inverse function of normal CDF
    alpha = (alpha < 0.001) ? 0.001 : alpha;
    return ppf[(int)(alpha * 1000) - 1];
}

double FValue(double alpha){//normal CDF
    alpha = (alpha < -3.1) ? -3.1 : alpha;
    alpha = (alpha > 3.1) ? 3.1 : alpha;
    int ind = (int)((alpha + 3.1) * 1000);
    return cdf[ind];
}

struct path{
    double first, second;//mean and variance
    int length = 1;//path length
    int ubmax = -1, lbmin = -1;//upper bound maximizer and lower bound minimizer
    path(double _first, double _second){
        first = _first;
        second = _second;
    }
};

struct Pe{
    int first;//u in X(v)
    vector<path> second;//Pe sets from v to u
    Pe(int _first, vector<path> &_second){
        first = _first;
        second = _second;
    }
};

struct node{
    int level;
    int ranks;
    int parent;
    vector<int> children;
    vector<Pe> X;
};
node T[MAX_V];

struct label{
    int first;//X(v)'s ancestor X(u)
    vector<path> second;//non-dominated path set from v to u
    double varmax = -DBL_MAX, varmin = DBL_MAX;
    label(int _first, vector<path> &_second){
        first = _first;
        second = _second;
        for (int i = 0; i < second.size();i++){
            varmax = max(second[i].second, varmax);
            varmin = min(second[i].second, varmin);
            double ub = -DBL_MAX, lb = DBL_MAX;
            for (int j = 0; j < second.size(); j++){
                if (i == j)
                    continue;
                double phiInter = FValue((second[i].first - second[j].first) / (sqrt(second[j].second) - sqrt(second[i].second)));
                if (j < i){
                    if (phiInter > ub){
                        ub = phiInter;
                        second[i].ubmax = j;
                    }
                }
                else if (j > i){
                    if (phiInter < lb){
                        lb = phiInter;
                        second[i].lbmin = j;
                    }
                }
            }
        }
    }
};
vector<label> L[MAX_V];
vector<int> anc[MAX_V];
int root = -1;
unordered_map<int,vector<path> > adj[MAX_V];
unordered_map<int,int> adjo[MAX_V];
vector<int> order;
bool flag[MAX_V];
bool cmprank(const Pe &a, const Pe &b){
    return T[a.first].ranks > T[b.first].ranks;
}

bool cmppath(const path &a, const path &b){
    return a.first < b.first;
}

vector<int> ordergen;
int del[MAX_V];//deleted neighbors
double fetchPriority(int v){
    return 1000 * adjo[v].size() + del[v];
}

double fetchTrav(double mean, double var, double alpha){
    return mean + FinvValue(alpha) * sqrt(var);
}

struct edge{
    int from, to;
    double mean, variance;
    edge(int a,int b,double c,double d){
        from = a, to = b, mean = c, variance = d;
    }
};
vector<edge> alledges;

void genorder(string filename, bool writeflag){
    priority_queue<II, vector<II>, greater<II> > degque;
    for (int i = 0; i < N; i++)
        degque.push(II(fetchPriority(i), i));
    int iter = -1, totnewedge = 0;
    while(!degque.empty()){
        II ii = degque.top();
        degque.pop();
        int v = ii.second;
        if(flag[v])
            continue;
        double prio = fetchPriority(v);
        if (prio > degque.top().first){
            degque.push(II(prio,v));
            continue;
        }
        iter += 1;
        //if(iter%100000==0){
        //    printf("----------%d %d--------\n", iter, totnewedge);
        //}
        flag[v] = 1;
        ordergen.push_back(v);
        T[v].ranks = iter;
        unordered_map<int, int>::iterator it;
        vector<int> nei;
        for (it = adjo[v].begin(); it !=adjo[v].end(); it++)
            if(!flag[it->first])
                nei.push_back(it->first);
        int lenX = nei.size();
        for (int j = 0; j < lenX; j++){
            int u = nei[j];
            for (int k = j + 1; k < lenX; k++){
                int w = nei[k];
                if(adjo[u].count(w)==0){
                    adjo[u][w] = 1;
                    adjo[w][u] = 1;
                    totnewedge += 1;
                }
            }
            //adjo[u].erase(v);
            del[u]++;
        }
    }
    if(writeflag){
        FILE *fp_order = fopen(filename.c_str(), "w");
        for (int i = 0; i < N;i++){
            fprintf(fp_order, "%d\n", T[i].ranks);
        }
        fclose(fp_order);
    }
}

void refine(vector<path> &P1,vector<path> &P2, vector<path> &res){
    if(P1.size()==0)
        for (int i = 0; i < P2.size();i++)
            res.push_back(P2[i]);
    if(P2.size()==0)
        for (int i = 0; i < P1.size();i++)
            res.push_back(P1[i]);
    for (int i = 0; i < P1.size();i++)
        for (int j = 0; j < P2.size();j++){
            res.push_back(path(P1[i].first + P2[j].first, P1[i].second + P2[j].second));
        }
    sort(res.begin(), res.end(), cmppath);
    double prev = DBL_MAX;
    vector<path> tmp;
    for (int i = 0; i < res.size();i++){
        if (res[i].first + sqrt(res[i].second) * 3.1 < prev){
            tmp.push_back(res[i]);
            prev = res[i].first + sqrt(res[i].second) * 3.1;
        }
    }
    res = tmp;
    //for (int i = 0; i < res.size();i++)
    //    printf("|%f %f|\n", res[i].first, res[i].second);
}

int descnt[MAX_V];
size_t avgPesize;
void treedec(){
    for (int i = 0; i < N; i++){
        int v = ordergen[i];
        //printf("|%d %d %d|\n", ii.first, adj[v].size(),v+1);
        if(i%100000==0)
            printf("%d\n", i);
        unordered_map<int, vector<path>>::iterator it;  
        for (it = adj[v].begin(); it !=adj[v].end(); it++)
            T[v].X.push_back(Pe(it->first, it->second));
        int lenX = T[v].X.size();
        /*for (int j = 0; j < lenX;j++){
            printf("%d ", T[v].X[j].first+1);
        }*/
        for (int j = 0; j < lenX; j++){
            Pe nu = T[v].X[j];
            int u = nu.first;
            for (int k = j + 1; k < lenX; k++){
                Pe nw = T[v].X[k];
                int w = nw.first;
                vector<path> emp;
                if(T[u].ranks<T[w].ranks){
                    if(adj[u].count(w)==0)
                        adj[u][w] = emp;
                    refine(nu.second, nw.second, adj[u][w]);
                }
                else{
                    if(adj[w].count(u)==0)
                        adj[w][u] = emp;
                    refine(nu.second, nw.second, adj[w][u]);
                }
            }
            avgPesize = max(avgPesize, nu.second.size());
        }
    }
    for (int i = 0; i < ordergen.size();i++){
        int v = ordergen[i];
        sort(T[v].X.begin(), T[v].X.end(), cmprank);
        int lenx = T[v].X.size();
        if (lenx != 0)
            T[v].parent = T[v].X[lenx - 1].first;
        else
            T[v].parent = MAX_V;
        //printf("%d %d\n", v + 1, T[v].parent + 1);
        vector<path> tmpd;
        T[v].X.push_back(Pe(v, tmpd));
        treewidth = max(treewidth, lenx + 1);
        if (T[v].parent == MAX_V){
            root = v;
            break;
        }
        T[T[v].parent].children.push_back(v);
        descnt[v]++;
        descnt[T[v].parent] += descnt[v];
    }
    cout << "Tree Width " << treewidth << endl;
}

vector<int> ancarray[MAX_V];
long long maxlabelsize, avglabelsize;
queue<int> bfs, bfssave;
void generateLabel4v(int v){
    vector<int> anc;
    int iu = v;
    //printf("v%d:", v+1);
    while (T[iu].parent != MAX_V){
        anc.push_back(T[iu].parent);
        iu = T[iu].parent;
    }
    int lenanc = anc.size();
    treeavgheight += lenanc;
    treeheight = max(treeheight, lenanc + 1);
    for (int i = 0; i < lenanc;i++){
        int u = anc[anc.size() - 1 - i];
        //printf("u%d ", u + 1);
        int lenx = T[v].X.size();
        vector<path> res;
        for (int j = 0; j < lenx;j++){
            int w = T[v].X[j].first;
            if (w == v)
                continue;
            //printf("w%d ", w + 1);
            if (T[w].ranks <= T[u].ranks){
                refine(T[v].X[j].second, L[w][i].second, res);
                if (w == u)
                    ancarray[v].push_back(i);
            }
            else{//w>u, w has been in j-th ancarray //ancarray[v][j] not used because need sorting ancarray
                refine(T[v].X[j].second, L[u][ancarray[v][j]].second, res);
            }
        }
        //for (int j = 0; j < res.size();j++)
        //    printf("res(%f,%f)", res[j].first, res[j].second);
        //vector<DD> tmp;
        //refine(tmp, tmp, res, sres);
        //cout << endl;
        //for (int j = 0; j < sres.size();j++)
        //    printf("sres%d(%f,%f)", u+1, sres[j].first, sres[j].second);
        L[v].push_back(label(u, res));
        maxlabelsize = max(maxlabelsize, (long long)res.size());
        avglabelsize += (long long)res.size();
        //cout << endl;
    }
    vector<path> tmpv;
    L[v].push_back(label(v, tmpv));
    ancarray[v].push_back(anc.size());
    //for (int j = 0; j < ancarray[v].size();j++)
    //    printf("anc%d ", ancarray[v][j]);
    //cout << L[v].size() << endl;
}
void labeling(){
    bfs.push(root);
    int iter = 0;
    while(!bfs.empty()){
        int v= bfs.front();
        bfs.pop();
        //sort(T[v].X.begin(), T[v].X.end(), cmp);
        generateLabel4v(v);
        for (int i = 0; i < T[v].children.size();i++)
            bfs.push(T[v].children[i]);
        if(iter%100000==0)
            printf("%d %d\n", iter, treeheight);
        iter += 1;
    }
}

long long indexsize;
void save(string filename){
    filename += string("index");
    ofstream of;
    of.open(filename.c_str(), ios::binary);
    // FILE *fp_index=fopen("index.txt","w");
    // fprintf(fp_index, "%d ", N);
    of.write(reinterpret_cast<const char *>(&N), sizeof(int));
    bfssave.push(root);
    while(!bfssave.empty()){
        int v = bfssave.front();
        bfssave.pop();
        //printf("%d\n", v);
        int lenl = L[v].size(), nx = T[v].X.size();
        indexsize = indexsize + 4 + nx;
        //fprintf(fp_index, "%d %d %d %d%c", v, T[v].parent, nx, lenl,' ');
        of.write(reinterpret_cast<const char *>(&v), sizeof(int));
        of.write(reinterpret_cast<const char *>(&T[v].parent), sizeof(int));
        of.write(reinterpret_cast<const char *>(&nx), sizeof(int));
        of.write(reinterpret_cast<const char *>(&lenl), sizeof(int));
        for (int i = 0; i < nx; i++){
            //fprintf(fp_index, "%d%c", T[v].X[i].first, (i == nx - 1) ? ' ' : ' ');
            of.write(reinterpret_cast<const char *>(&T[v].X[i].first), sizeof(int));
        }
        for (int i = 0; i < lenl; i++){
            int lend = L[v][i].second.size();
            indexsize = indexsize + 2 + lend * 2;
            //fprintf(fp_index, "%d %d ", L[v][i].first, lend);
            of.write(reinterpret_cast<const char *>(&L[v][i].first), sizeof(int));
            of.write(reinterpret_cast<const char *>(&lend), sizeof(int));
            for (int j = 0; j < lend;j++){
                //fprintf(fp_index, "%d %d ", L[v][i].second[j].first*10, L[v][i].second[j].second*10);
                int p = L[v][i].second[j].first, q = L[v][i].second[j].second;
                of.write(reinterpret_cast<const char *>(&p), sizeof(int));
                of.write(reinterpret_cast<const char *>(&q), sizeof(int));
            }
        }
        for (int i = 0; i < T[v].children.size();i++)
            bfssave.push(T[v].children[i]);
    }
    //fclose(fp_index);
    of.close();
}

void concatQuery(int s, int t, int ind, double alpha){
    vector<path> P1 = L[s][ind].second, P2 = L[t][ind].second;
    vector<bool> flag1(P1.size()), flag2(P2.size());
    double varmax1 = L[s][ind].varmax, varmin1 = L[s][ind].varmin;
    double varmax2 = L[t][ind].varmax, varmin2 = L[t][ind].varmin;
    for (int i = 0; i < P1.size();i++){
        int ubm = P1[i].ubmax, lbm = P1[i].lbmin, flag = 0;
        if(ubm!=-1){
            double ub = FValue((P1[i].first - P1[ubm].first) / (sqrt(P1[ubm].second + varmin2) - sqrt(P1[i].second + varmin2)));
            if (alpha < ub)
                flag1[i] = 1;
        }
        if (!flag && lbm != -1){
            double lb = FValue((P1[i].first - P1[lbm].first) / (sqrt(P1[lbm].second + varmax2) - sqrt(P1[i].second + varmax2)));
            if (alpha > lb)
                flag1[i] = 1;
        }
    }
    for (int i = 0; i < P2.size();i++){
        int ubm = P2[i].ubmax, lbm = P2[i].lbmin, flag = 0;
        if(ubm!=-1){
            double ub = FValue((P2[i].first - P2[ubm].first) / (sqrt(P2[ubm].second + varmin1) - sqrt(P2[i].second + varmin1)));
            if (alpha < ub)
                flag2[i] = 1;
        }
        if (!flag && lbm != -1){
            double lb = FValue((P2[i].first - P2[lbm].first) / (sqrt(P2[lbm].second + varmax1) - sqrt(P2[i].second + varmax1)));
            if (alpha > lb)
                flag2[i] = 1;
        }
    }
    for (int i = 0; i < P1.size();i++){
        if(flag1[i])
            continue;
        for (int j = 0; j < P2.size();j++){
            if(flag2[j])
                continue;
            npathConcat++;
            optw = min(optw, fetchTrav(P1[i].first + P2[j].first, P1[i].second + P2[j].second, alpha));
        }
	}
}

void NRPQuery(int s, int t, double alpha){
    optw = DBL_MAX;
    if (s == t)
        return;
    s--;
    t--;
    vector<int> ancs,anct;
    int u1 = s, l = -1;
    while(u1!=MAX_V){
        ancs.push_back(u1);
        u1 = T[u1].parent;
    }
    u1 = t;
    while(u1!=MAX_V){
        anct.push_back(u1);
        u1 = T[u1].parent;
    }
    int i = ancs.size() - 1, j = anct.size() - 1, k = -1;
    unordered_map<int, int> inds, indt;
    while (i != -1 && j != -1){
        if (ancs[i] == anct[j]){
            inds[ancs[i]] = i;
            indt[anct[j]] = j;
            i--;
            j--;
            k++;
        }
        else
            break;
    }
    if(i==-1)
        l = ancs[0];
    else if(j==-1)
        l = anct[0];
    else
        l = ancs[i + 1];
    int ind;
    if (l == s){
        int lenp = L[t][ancs.size() - 1].second.size();
        for (int i = 0; i < lenp;i++){
            double trav=fetchTrav(L[t][ancs.size() - 1].second[i].first, L[t][ancs.size() - 1].second[i].second, alpha);
            optw = min(optw, trav);
        }
    }
    else if (l == t){
        int lenp = L[s][anct.size() - 1].second.size();
        for (int i = 0; i < lenp;i++){
            double trav=fetchTrav(L[s][anct.size() - 1].second[i].first, L[s][anct.size() - 1].second[i].second, alpha);
            optw = min(optw, trav);
        }
    }
    else{
        int cs = ancs[i], ct = anct[j];
        l = (ancarray[cs].size() < ancarray[ct].size()) ? cs : ct;
        //printf("lca%d\n", l + 1);
        hopsize += ancarray[l].size() - 1;
        for (int i = 0; i < ancarray[l].size() - 1; i++){
            ind = ancarray[l][i];
            //printf("|%d %d %d|\n", ind, L[s][ind].second.size(), L[t][ind].second.size());
            concatQuery(s, t, ind, alpha);
        }
    }
    //printf("%f\n", optw);
}

int main(int argc , char * argv[]){
    string sfile, sq;
    FILE *fp_query, *fp_networkvar, *fp_networkmu;
    if (argc > 1)
        sfile = string(argv[1]);
    else
        sfile = string("NY");
    if (argc > 2)
        sq = string(argv[2]);
    else
        sq = string("pool1");
    string prefix = string("../data/") + sfile + string("/");
    string smu = prefix + string("USA-road-t.") + sfile + (".gr");
    string svar = prefix + string("USA-road-var.") + sfile + (".gr");
    //test or not
    if(0){
        sq = string("test");
        smu = string("../data/") + sfile + string("/t.") + sfile + (".gr");
        svar = string("../data/") + sfile + string("/var.") + sfile + (".gr");
    }
    fp_networkmu = fopen(smu.c_str(), "r");
    fp_networkvar = fopen(svar.c_str(), "r");
    char ch, buffer[100];
    int u, v;
    double w, c;
    //travel time
    for (int i = 0; i < 4; i++){
        fgets(buffer, 90, fp_networkvar);
        fgets(buffer, 90, fp_networkmu);
    }
    for (int i = 0; i < 4; i++){
        fgetc(fp_networkvar);
        fgetc(fp_networkmu);
    }
    fscanf(fp_networkvar, "%d%d", &N, &M);
    fscanf(fp_networkmu, "%d%d", &N, &M);
    for (int i = 0; i < 3; i++){
        fgets(buffer, 90, fp_networkvar);
        fgets(buffer, 90, fp_networkmu);
    }
    for (int i = 0; i < M; i++) {
        fscanf(fp_networkvar, "%c%d%d%lf", &ch, &u, &v, &w);
        fgets(buffer, 90, fp_networkvar);
        fscanf(fp_networkmu, "%c%d%d%lf", &ch, &u, &v, &c);
        fgets(buffer, 90, fp_networkmu);
        u--;
        v--;        
        if(i%2==0){
            adjo[u][v]=1;
            adjo[v][u]=1;
            alledges.push_back(edge(u, v, c, w));
        }
    }

    std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT;
    //generate order for vertices
    string ordername = string("../data/") + sfile + string("/") + string("order.txt");
    if(0){
        t1=std::chrono::high_resolution_clock::now();
        genorder(ordername, 0);
        t2 = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        runT= time_span.count();
        cout<<"Order Generation Time "<<runT<<endl;
    }
    else{
        ordergen.assign(N, -1);
        FILE *fpord = fopen(ordername.c_str(), "r");
        for (int i = 0; i < N; i++){
            fscanf(fpord, "%d", &T[i].ranks);
            ordergen[T[i].ranks] = i;
        }
    }
    
    //distribute edges
    for (int i = 0; i < alledges.size();i++){
        int f = alledges[i].from, t = alledges[i].to;
        if(T[f].ranks<T[t].ranks)
            adj[f][t].push_back(path(alledges[i].mean, alledges[i].variance));
        else
            adj[t][f].push_back(path(alledges[i].mean, alledges[i].variance));
    }

    t1=std::chrono::high_resolution_clock::now();
    treedec();
    t2=std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"Tree Decomposition Time "<<runT<<endl;
    //return 0;

    t1=std::chrono::high_resolution_clock::now();
    labeling();
    t2=std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"Labeling Time "<<runT<<endl;
    cout<<"Tree Height "<<treeheight<<endl;
    cout<<"Tree Avg. Height "<<(double)treeavgheight/N<<endl;
    cout << "Max. Label Size " << maxlabelsize << endl;
    cout << "Avg. Label Size " << (double)avglabelsize / treeavgheight << endl;    

    
    t1=std::chrono::high_resolution_clock::now();
    //FILE *fp = fopen("index", "w");
    //save(string("../data/") + sfile + string("/"));
    t2=std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"Saving Time "<<runT<<endl;
    cout << "Index Size " << (double)indexsize * 4 / 1000000 << "MB" << endl;

    freopen((prefix + string("Results_ind")).c_str(), "w", stdout);
    for (int i = 0; i < 5;i++){
        hopsize = npathConcat = 0;
        string s3 = string("../data/") + sfile + string("/") + string("q") + to_string(i+1);
        fp_query = fopen(s3.c_str(), "r");
        vector<pair<II, double>> queryset;
        int qs, qt;
        double qC;
        while (~fscanf(fp_query, "%d%d%lf", &qs, &qt, &qC)){
            queryset.push_back(make_pair(II(qs, qt), qC));
        }
        FILE *fp_ans1 = fopen((s3+string("RSPans_ind")).c_str(), "w");
        t1=std::chrono::high_resolution_clock::now();
        for(int i=0;i<queryset.size();i++){
            NRPQuery(queryset[i].first.first, queryset[i].first.second, queryset[i].second);
            fprintf(fp_ans1, "%f\n", optw);
        }
        fclose(fp_ans1);
        t2=std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        runT= time_span.count();
        cout<<"Query Time "<<runT<<endl;
        cout << "# of Hoplinks " << hopsize <<endl;
        cout << "# of Path Concatenations " << npathConcat <<endl;
    }

    for (int i = 0; i < 5;i++){
        hopsize = npathConcat = 0;
        string s3 = string("../data/") + sfile + string("/") + string("r") + to_string(i+1);
        fp_query = fopen(s3.c_str(), "r");
        vector<pair<II, double>> queryset;
        int qs, qt;
        double qC;
        while (~fscanf(fp_query, "%d%d%lf", &qs, &qt, &qC)){
            queryset.push_back(make_pair(II(qs, qt), qC));
        }
        FILE *fp_ans1 = fopen((s3+string("RSPans_ind")).c_str(), "w");
        t1=std::chrono::high_resolution_clock::now();
        for(int i=0;i<queryset.size();i++){
            NRPQuery(queryset[i].first.first, queryset[i].first.second, queryset[i].second);
            fprintf(fp_ans1, "%f\n", optw);
        }
        fclose(fp_ans1);
        t2=std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        runT= time_span.count();
        cout<<"Query Time "<<runT<<endl;
        cout << "# of Hoplinks " << hopsize <<endl;
        cout << "# of Concatenations " << npathConcat <<endl;
    }
}