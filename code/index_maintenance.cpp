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
const int MAX_N = 435676, MAX_M = 1057266;
int neiK = 5;
int N, M;
long long npathConcat, hopsize, p12size;
double optw;
int treeheight = 0, treewidth = 0, treeavgheight = 0;
unordered_map<int, double> cov[MAX_M];
bool corrFlag[MAX_N];

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
    int f, t;
    vector<int> from, to;// head from the vertex of lower rank ("from"); tail from "to"
    int ubmax = -1, lbmin = -1;//upper bound maximizer and lower bound minimizer    
    path(int _f, int _t, vector<int> &_from, vector<int> &_to, double _first, double _second){
        f = _f;
        t = _t;
        from = _from;
        to = _to;
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
node T[MAX_N];

struct label{
    int first;//u: X(v)'s ancestor X(u)
    vector<path> second;//non-dominated path set from v to u
    double varmax = -DBL_MAX, varmin = DBL_MAX;
    label(int _first, vector<path> &_second){
        first = _first;
        second = _second;
        for (int i = 0; i < second.size();i++){
            varmax = max(second[i].second, varmax);
            if(corrFlag[first]==0){//only u will be the hoplink
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
    }
};
vector<label> L[MAX_N];
int root = -1;
unordered_map<int, vector<path>> adj[MAX_N]; // during tree decomposition
unordered_map<int, int> adjo[MAX_N];
unordered_map<int, vector<int>> C[MAX_N];
vector<int> order;
bool cmprank(const Pe &a, const Pe &b){
    return T[a.first].ranks > T[b.first].ranks;
}

bool cmppath(const path &a, const path &b){
    return a.first < b.first;
}

vector<int> ordergen;
int del[MAX_N];//deleted neighbors
double fetchPriority(int v){
    return 1000 * adjo[v].size() + del[v];
}

double fetchTrav(double mean, double var, double alpha){
    return mean + FinvValue(alpha) * sqrt(var);
}

struct edge{
    int from, to, id;
    double mean, variance;
    edge(int a,int b,double c,double d, int _id){
        from = a, to = b, mean = c, variance = d, id = _id;
    }
};
vector<edge> alledges;

bool flagOrder[MAX_N];
void genorder(string filename, bool writeflag){
    priority_queue<II, vector<II>, greater<II> > degque;
    for (int i = 0; i < N; i++)
        degque.push(II(fetchPriority(i), i));
    int iter = -1, totnewedge = 0;
    while(!degque.empty()){
        II ii = degque.top();
        degque.pop();
        int v = ii.second;
        if(flagOrder[v])
            continue;
        double prio = fetchPriority(v);
        if (prio > degque.top().first){
            degque.push(II(prio,v));
            continue;
        }
        iter += 1;
        flagOrder[v] = 1;
        ordergen.push_back(v);
        T[v].ranks = iter;
        unordered_map<int, int>::iterator it;
        vector<int> nei;
        for (it = adjo[v].begin(); it !=adjo[v].end(); it++)
            if(!flagOrder[it->first])
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
path combine2Paths(path &p1, path &p2){
    double variance = p1.second + p2.second;
    int p1f = p1.f, p2f = p2.f, p1t = p1.t, p2t = p2.t, f, t;
    int piv;
    vector<int> *newf, *newt, *cov1, *cov2;
    if(p1f==p2f){//tree decomposition
        piv = p1f;
        if(T[p1t].ranks<T[p2t].ranks)
            f = p1t, t = p2t, newf = &p1.to, newt = &p2.to, cov1 = &p2.from, cov2 = &p1.from;
        else
            f = p2t, t = p1t, newf = &p2.to, newt = &p1.to, cov1 = &p1.from, cov2 = &p2.from;
    }
    else if(p1f==p2t){
        piv = p1f;
        if(T[p1t].ranks<T[p2f].ranks)
            f = p1t, t = p2f, newf = &p1.to, newt = &p2.from, cov1 = &p2.to, cov2 = &p1.from;
        else
            f = p2f, t = p1t, newf = &p2.from, newt = &p1.to, cov1 = &p1.from, cov2 = &p2.to;
    }
    else if(p1t==p2f){//label v w u where T[w].ranks<T[u].ranks
        piv = p1t;
        if(T[p1f].ranks<T[p2t].ranks)
            f = p1f, t = p2t, newf = &p1.from, newt = &p2.to, cov1 = &p2.from, cov2 = &p1.to;
        else
            f = p2t, t = p1f, newf = &p2.to, newt = &p1.from, cov1 = &p1.to, cov2 = &p2.from;
    }
    else if(p1t==p2t){//label v w u where T[u].ranks<T[w].ranks
        piv = p1t;
        if(T[p1f].ranks<T[p2f].ranks)
            f = p1f, t = p2f, newf = &p1.from, newt = &p2.from, cov1 = &p2.to, cov2 = &p1.to;
        else
            f = p2f, t = p1f, newf = &p2.from, newt = &p1.from, cov1 = &p1.to, cov2 = &p2.to;
    }
    else{
        printf("%d %d %d %d improper path combination\n", p1f, p1t, p2f, p2t);
    }
    if(corrFlag[piv]){
        for (int i = 0; i < cov1->size();i++)
            for (int j = 0; j < cov2->size();j++){
                variance += cov[(*cov1)[i]][(*cov2)[j]] * 2;
            }
    }
    if (variance < 0)
        printf("variance of combined path less than 0\n");
    path ret = path(f, t, *newf, *newt, p1.first + p2.first, variance);
    size_t rm1 = neiK - newf->size(), rm2 = neiK - newt->size();
    for (int i = 0; i < min(cov1->size(),rm1); i++)
        ret.from.push_back((*cov1)[i]);
    for (int i = 0; i < min(cov2->size(),rm2); i++)
        ret.to.push_back((*cov2)[i]);    
    ret.length = p1.length + p2.length;
    return ret;
}

bool flagedges[MAX_M] = {0};
bool checkdom(path &p1, path &p2){
    bool domflag = 1;
    if(corrFlag[p1.f]){
        for (int i = 0; i < p1.from.size();i++)
            flagedges[p1.from[i]] = 1;
        for (int i = 0; i < p2.from.size();i++)
            flagedges[p2.from[i]] = 1;
        unordered_map<int, int>::iterator it;
        for (it = adjo[p1.f].begin(); it != adjo[p1.f].end(); it++){
            int eid = it->second;
            if(flagedges[eid])
                continue;
            double var1 = p1.second + alledges[eid].variance;
            double var2 = p2.second + alledges[eid].variance;
            for (int i = 0; i < p1.from.size();i++)
                var1 += 2*cov[eid][p1.from[i]];
            for (int i = 0; i < p2.from.size();i++)
                var2 += 2*cov[eid][p2.from[i]];
            if(var1<0||var2<0)
                printf("var12<0\n");
            if(p1.first+3.1*sqrt(var1)>p2.first+3.1*sqrt(var2))
                domflag = 0;
        }
        for (int i = 0; i < p1.from.size();i++)
            flagedges[p1.from[i]] = 0;
        for (int i = 0; i < p2.from.size();i++)
            flagedges[p2.from[i]] = 0;
    }
    if(corrFlag[p1.t]){
        for (int i = 0; i < p1.to.size();i++)
            flagedges[p1.to[i]] = 1;
        for (int i = 0; i < p2.to.size();i++)
            flagedges[p2.to[i]] = 1;
        unordered_map<int, int>::iterator it;
        for (it = adjo[p1.t].begin(); it != adjo[p1.t].end(); it++){
            int eid = it->second;
            if(flagedges[eid])
                continue;
            double var1 = p1.second + alledges[eid].variance;
            double var2 = p2.second + alledges[eid].variance;
            for (int i = 0; i < p1.to.size();i++)
                var1 += 2*cov[eid][p1.to[i]];
            for (int i = 0; i < p2.to.size();i++)
                var2 += 2*cov[eid][p2.to[i]];
            if(var1<0||var2<0)
                printf("var12<0\n");
            if(p1.first+3.1*sqrt(var1)>p2.first+3.1*sqrt(var2))
                domflag = 0;
        }
        for (int i = 0; i < p1.to.size();i++)
            flagedges[p1.to[i]] = 0;
        for (int i = 0; i < p2.to.size();i++)
            flagedges[p2.to[i]] = 0;
    }
    return domflag;
}

void refine(vector<path> &P1,vector<path> &P2, vector<path> &res){
    if(P1.size()==0)
        for (int i = 0; i < P2.size();i++)
            res.push_back(P2[i]);
    if(P2.size()==0)
        for (int i = 0; i < P1.size();i++)
            res.push_back(P1[i]);
    double prev = DBL_MAX;
    if(P1.size()>P2.size()){
        vector<bool> P1F(P1.size(), 0);
        for (int i = 0; i < P1.size();i++){
            if (P1[i].first + sqrt(P1[i].second) * 3.1 < prev){
                P1F[i] = 1;
                prev = P1[i].first + sqrt(P1[i].second) * 3.1;
            }
        }
        for (int j = 0; j < P2.size();j++){
            if(P2[j].length>=neiK)
                for (int i = 0; i < P1.size();i++){
                    if(P1F[i])
                        res.push_back(combine2Paths(P1[i], P2[j]));
                }
            else
                for (int i = 0; i < P1.size();i++)
                    res.push_back(combine2Paths(P1[i], P2[j]));
        }
    }
    else{
        vector<bool> P2F(P2.size(), 0);
        for (int i = 0; i < P2.size();i++){
            if (P2[i].first + sqrt(P2[i].second) * 3.1 < prev){
                P2F[i] = 1;
                prev = P2[i].first + sqrt(P2[i].second) * 3.1;
            }
        }
        for (int i = 0; i < P1.size();i++){
            if(P1[i].length>=neiK)
                for (int j = 0; j < P2.size();j++){
                    if(P2F[j])
                        res.push_back(combine2Paths(P1[i], P2[j]));
                }
            else
                for (int j = 0; j < P2.size();j++)
                    res.push_back(combine2Paths(P1[i], P2[j]));
        }
    }    
    /*
    for (int i = 0; i < P1.size();i++)
        for (int j = 0; j < P2.size();j++){
            res.push_back(combine2Paths(P1[i], P2[j]));
        }
    */
    sort(res.begin(), res.end(), cmppath);
    prev = DBL_MAX;
    vector<path> tmp;
    int from = res[0].f, to = res[0].t;
    if(!corrFlag[from]&&!corrFlag[to]){
        for (int i = 0; i < res.size();i++){
            if (res[i].first + sqrt(res[i].second) * 3.1 < prev){
                tmp.push_back(res[i]);
                prev = res[i].first + sqrt(res[i].second) * 3.1;
            }
        }
    }
    else{
        tmp.push_back(res[0]);
        for (int i = 1; i < res.size();i++){
            if (!checkdom(tmp[tmp.size() - 1], res[i]))
                tmp.push_back(res[i]);
        }
    }
    res = tmp;
    //for (int i = 0; i < res.size();i++)
    //    printf("|%f %f|\n", res[i].first, res[i].second);
}

int descnt[MAX_N];
size_t avgPesize;
void treedec(){
    for (int i = 0; i < N; i++){
        int v = ordergen[i];
        if(i%100000==0)
            printf("%d %d\n", i, avgPesize);
        if(i>(N/100000)*100000&&i%10000==0)
            printf("%d %d\n", i, avgPesize);
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
                    C[u][w].push_back(v);
                }
                else{
                    if(adj[w].count(u)==0)
                        adj[w][u] = emp;
                    refine(nu.second, nw.second, adj[w][u]);
                    C[w][u].push_back(v);
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
            T[v].parent = MAX_N;
        vector<path> tmpd;
        T[v].X.push_back(Pe(v, tmpd));
        treewidth = max(treewidth, lenx + 1);
        if (T[v].parent == MAX_N){
            root = v;
            break;
        }
        T[T[v].parent].children.push_back(v);
        descnt[v]++;
        descnt[T[v].parent] += descnt[v];
    }
}

vector<int> ancarray[MAX_N];
long long maxlabelsize, avglabelsize;
queue<int> bfs, bfssave;
void generateLabel4v(int v){
    vector<int> anc;
    int iu = v;//printf("v%d:", v+1);
    while (T[iu].parent != MAX_N){
        anc.push_back(T[iu].parent);
        iu = T[iu].parent;
    }
    int lenanc = anc.size();
    treeavgheight += lenanc;
    treeheight = max(treeheight, lenanc + 1);
    for (int i = 0; i < lenanc;i++){
        int u = anc[anc.size() - 1 - i];//printf("u%d ", u + 1);
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
    filename += string("NRPindex");
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
            indexsize = indexsize + 2 + lend * (3 + 2 * 2);
            //fprintf(fp_index, "%d %d ", L[v][i].first, lend);
            of.write(reinterpret_cast<const char *>(&L[v][i].first), sizeof(int));
            of.write(reinterpret_cast<const char *>(&lend), sizeof(int));
            for (int j = 0; j < lend;j++){
                //fprintf(fp_index, "%d %d ", L[v][i].second[j].first*10, L[v][i].second[j].second*10);
                int fr = L[v][i].second[j].f, to = L[v][i].second[j].t, l = L[v][i].second[j].length;
                of.write(reinterpret_cast<const char *>(&fr), sizeof(int));
                of.write(reinterpret_cast<const char *>(&to), sizeof(int));
                of.write(reinterpret_cast<const char *>(&l), sizeof(int));
                int p = L[v][i].second[j].first, q = L[v][i].second[j].second;
                of.write(reinterpret_cast<const char *>(&p), sizeof(double));
                of.write(reinterpret_cast<const char *>(&q), sizeof(double));
                for (int k = 0; k < L[v][i].second[j].from.size();k++)
                    of.write(reinterpret_cast<const char *>(&L[v][i].second[j].from[k]), sizeof(int));
                for (int k = 0; k < L[v][i].second[j].to.size();k++)
                    of.write(reinterpret_cast<const char *>(&L[v][i].second[j].to[k]), sizeof(int));
                indexsize += L[v][i].second[j].from.size() + L[v][i].second[j].to.size();
            }
        }
        for (int i = 0; i < T[v].children.size();i++)
            bfssave.push(T[v].children[i]);
    }
    //fclose(fp_index);
    of.close();
}

bool updateflagv[MAX_N];
unordered_map<int,int> updateE[MAX_N];
void updatePuw(int u, int w){
    int lenC = C[u][w].size();
    vector<path> emp;
    int uwpos = lower_bound(T[u].X.begin(), T[u].X.end(), Pe(w, emp), cmprank) - T[u].X.begin();
    vector<DD> original;
    int lenux = T[u].X[uwpos].second.size();
    for (int i = 0; i < lenux; i++)
        original.push_back(DD(T[u].X[uwpos].second[i].first, T[u].X[uwpos].second[i].second));
    T[u].X[uwpos].second = emp;
    if(adjo[u].count(w)){
        int uweid = adjo[u][w];
        vector<int> p;
        if (neiK > 0)
            p.push_back(uweid);
        T[u].X[uwpos].second.push_back(path(u, w, p, p, alledges[uweid].mean, alledges[uweid].variance));
    }
    for (int i = 0; i < lenC; i++){
        int v = C[u][w][i];
        int vupos = lower_bound(T[v].X.begin(), T[v].X.end(), Pe(u,emp), cmprank) - T[v].X.begin();
        int vwpos = lower_bound(T[v].X.begin(), T[v].X.end(), Pe(w,emp), cmprank) - T[v].X.begin();
        refine(T[v].X[vupos].second, T[v].X[vwpos].second, T[u].X[uwpos].second);
    }

    bool uwchanges = 0;
    if (lenux != T[u].X[uwpos].second.size())
        uwchanges = 1;
    else{
        for (int i = 0; i < lenux; i++){
            if (abs(T[u].X[uwpos].second[i].first - original[i].first) > 1e-9)
                uwchanges = 1;
            if (abs(T[u].X[uwpos].second[i].second - original[i].second) > 1e-9)
                uwchanges = 1;
            if(uwchanges)
                break;
        }
    }
    
    if (uwchanges == 0){
        updateflagv[u] = 1;
        return;
    }
    else{
        updateflagv[u] = 0;
    }    

    int lenx = T[u].X.size();
    for (int i = lenx - 1; i >= 0; i--){
        int v = T[u].X[i].first;
        if (v == u || v == w)
            continue;
        if(T[v].ranks<T[w].ranks){
            updateflagv[v] = 1;
            updateE[v][w] = 1;
        }
        else{
            updateflagv[w] = 1;
            updateE[w][v] = 1;
        }
    }
}

bool cmperank(const edge &a, const edge &b){
    return T[a.from].ranks < T[b.from].ranks;
}
bool cmpnrank(const int a, const int b){
    return T[a].ranks > T[b].ranks;
}

void updateLabel4v(int v){
    vector<int> anc;
    int iu = v;
    while (T[iu].parent != MAX_N){
        anc.push_back(T[iu].parent);
        iu = T[iu].parent;
    }
    int lenanc = anc.size();
    for (int i = 0; i < lenanc;i++){
        int u = anc[anc.size() - 1 - i];
        int lenx = T[v].X.size();
        vector<path> res;
        for (int j = 0; j < lenx;j++){
            int w = T[v].X[j].first;
            if (w == v)
                continue;
            if (T[w].ranks <= T[u].ranks){
                refine(T[v].X[j].second, L[w][i].second, res);
            }
            else{//w>u, w has been in j-th ancarray //ancarray[v][j] not used because need sorting ancarray
                refine(T[v].X[j].second, L[u][ancarray[v][j]].second, res);
            }
        }
        L[v][i] = label(u, res);
        maxlabelsize = max(maxlabelsize, (long long)res.size());
    }
}
void indexMaintenance(vector<edge> &newedges){
    for (int i = 0; i < newedges.size();i++){
        int f = newedges[i].from, t = newedges[i].to, tmp, eid = newedges[i].id;
        alledges[eid].mean = newedges[i].mean;
        alledges[eid].variance = newedges[i].variance;
        if(T[f].ranks>T[t].ranks){
            updateflagv[t] = 1;
            updateE[t][f] = 1;
        }
        else{
            updateflagv[f] = 1;
            updateE[f][t] = 1;
        }
    }

    std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT;
    t1=std::chrono::high_resolution_clock::now();    
    for (int i = 0; i < N; i++){
        int v = ordergen[i];
        if(updateflagv[v]){
            unordered_map<int, int>::iterator it;
            for (it = updateE[v].begin(); it != updateE[v].end(); it++)
                updatePuw(v, it->first);
        }
    }

    t2=std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"Pe Update Time "<<runT/1000<<endl;

    vector<int> vuroot;
    for (int i = 0; i < N;i++)
        if(updateflagv[i])
            vuroot.push_back(i);
    t1=std::chrono::high_resolution_clock::now();
    bool labelflag[MAX_N]={0};
    sort(vuroot.begin(), vuroot.end(), cmpnrank);
    printf("Labelupdate vurootsize %d %d\n", vuroot.size(), T[vuroot[0]].ranks);
    for (int i = 0; i < vuroot.size();i++){
        queue<int> vbfs;
        int uroot = vuroot[i];
        if(labelflag[uroot])
            continue;
        vbfs.push(uroot);
        while(!vbfs.empty()){
            int v = vbfs.front();
            vbfs.pop();
            if(labelflag[v])
                continue;
            updateLabel4v(v);
            labelflag[v] = 1;
            for (int i = 0; i < T[v].children.size();i++)
                vbfs.push(T[v].children[i]);
        }
    }
    t2=std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"Label Update Time "<<runT/1000<<endl;
}

int main(int argc , char * argv[]){
    string sfile, sq, st, sup;
    FILE *fp_query, *fp_networkvar, *fp_networkmu, *fp_networkcov, *fp_udpate;
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
        if(argv[3][0]=='K')
        neiK = atoi(st.substr(1).c_str());
    }
    else
        st = string("");
    string prefix = string("../data/") + sfile + string("/");
    string smu = prefix + string("USA-road-t.") + sfile + (".gr");
    string svar = prefix + sq+string("USA-road-var.") + sfile + (".gr");
    string scov = prefix + st+string("USA-road-cov.") + sfile + (".gr");
    if(0){//test or not
        sq = string("test");
        smu = string("../data/") + sfile + string("/t.") + sfile + (".gr");
        svar = string("../data/") + sfile + string("/var.") + sfile + (".gr");
        scov = string("../data/") + sfile + string("/cov.") + sfile + (".gr");
    }
    fp_networkmu = fopen(smu.c_str(), "r");
    fp_networkvar = fopen(svar.c_str(), "r");
    fp_networkcov = fopen(scov.c_str(), "r");
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
        fscanf(fp_networkmu, "%c%d%d%lf", &ch, &u, &v, &c);
        fgets(buffer, 90, fp_networkmu);
        fscanf(fp_networkvar, "%c%d%d%lf", &ch, &u, &v, &w);
        fgets(buffer, 90, fp_networkvar);
        u--;
        v--;        
        if(i%2==0){
            adjo[u][v] = i/2;
            adjo[v][u] = i/2;
            alledges.push_back(edge(u, v, c, w, i/2));
        }
    }
    while (fscanf(fp_networkcov, "%d%d%lf", &u, &v, &c) == 3){
        cov[u][v] = cov[v][u] = c;
        corrFlag[alledges[u].from] = corrFlag[alledges[u].to] = 1;
        corrFlag[alledges[v].from] = corrFlag[alledges[v].to] = 1;
    }
    printf("Cov Read\n");

    std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double runT;
    string ordername = string("../data/") + sfile + string("/") + string("order.txt");
    if(0){//generate order for vertices
        t1=std::chrono::high_resolution_clock::now();
        genorder(ordername, 0);
        t2 = std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        runT= time_span.count();
        cout<<"Order Generation Time "<<runT<<endl;
    }
    else{//get order from file
        ordergen.assign(N, -1);
        FILE *fpord = fopen(ordername.c_str(), "r");
        for (int i = 0; i < N; i++){
            fscanf(fpord, "%d", &T[i].ranks);
            ordergen[T[i].ranks] = i;
        }
    }
    
    //distribute edges to adj according to order
    for (int i = 0; i < alledges.size();i++){
        int f = alledges[i].from, t = alledges[i].to;
        vector<int> p;
        if(neiK>0)
            p.push_back(i);
        if(T[f].ranks<T[t].ranks)
            adj[f][t].push_back(path(f,t,p,p, alledges[i].mean, alledges[i].variance));
        else
            adj[t][f].push_back(path(f,t,p,p, alledges[i].mean, alledges[i].variance));
    }

    string setres;

    t1=std::chrono::high_resolution_clock::now();
    treedec();
    t2=std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"Tree Decomposition Time "<<runT<<endl;
    cout << "Tree Width " << treewidth << endl;
    setres += string("Tree Decomposition Time ") + to_string(runT)+ string("\n");

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
    setres += string("Labeling Time ") + to_string(runT) + string("\n");
    setres += string("Tree Width ") + to_string(treewidth) + string("\n");    
    setres += string("Tree Height ") + to_string(treeheight) + string("\n");   
    setres += string("Tree Avg. Height ") + to_string((double)treeavgheight/N) + string("\n");
    setres += string("Max. Label Size ") + to_string(maxlabelsize) + string("\n");
    setres += string("Avg. Label Size ") + to_string((double)avglabelsize / treeavgheight) + string("\n");   

    
    t1=std::chrono::high_resolution_clock::now();
    //save(string("../data/") + sfile + string("/"));
    t2=std::chrono::high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
	runT= time_span.count();
	cout<<"Saving Time "<<runT<<endl;
    cout << "Index Size " << (double)indexsize * 4 / 1000000 << "MB" << endl;
    setres += string("Saving Time ") + to_string(runT) + string("\n");
    setres += string("Index Size ") + to_string((double)indexsize * 4 / 1000000) + string("MB\n");

    freopen((prefix + string("updateResults")).c_str(), "a", stdout);
    if (argc > 4){
        cout << sfile << endl;
        fp_udpate = fopen((prefix + string(argv[4])).c_str(), "r");
        vector<edge> newedges;
        while (fscanf(fp_udpate, "%d%lf%lf", &u, &w, &c) == 3){
            int id = u;
            newedges.push_back(edge(alledges[id].from, alledges[id].to, alledges[id].mean * w, alledges[id].variance * c, id));
        }
        t1=std::chrono::high_resolution_clock::now();
        indexMaintenance(newedges);
        t2=std::chrono::high_resolution_clock::now();
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
        runT= time_span.count();
        cout<<"Update Time "<<runT/newedges.size()<<endl;

        int storage = 0;
        for (int i = 0; i < N;i++){
            storage += 1 + C[i].size();
            unordered_map<int, vector<int>>::iterator it;
            for (it = C[i].begin(); it != C[i].end();it++){
                storage += it->second.size();
            }
        }
        cout << "Extra Storage " << (double)storage * 4 / 1000000 << "MB" << endl;
    }

    return 0;
}