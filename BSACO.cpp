#include<iostream>
#include<math.h>
#include<set>

using namespace std;

const int N = 50;
const int M = 10;

long long seed = 131;

struct Point{
    int x, y;
};
struct PointSet{
    int siz;
    int brk1, brk2;
    double sum;
    vector<int>p;
    vector<int>subpth;
};

Point tp[N];
set<int>lst;
vector<PointSet>ps;

double dist(int a, int b){
    return sqrt((tp[a].x - tp[b].x) * (tp[a].x - tp[b].x) + (tp[a].y - tp[b].y) * (tp[a].y - tp[b].y));
}

namespace ACO{
    const int ant_num = 10;
    const double alpha = 1;
    const double beta = 5;
    const double rho = 0.1;
    const double delta = 1e-3;
    const double Q = 100;
    const int iter_max = 200;
    int city_num;
    double Tau[M][M], dis[M][M], Eta[M][M], DeltaTau[M][M];
    double Relation[M][M][M];
    double pathlen[N];
    int AntPath[N][M];
    int c[M];
    double P[M];

    double cossim(Point i, Point j, Point k){
        double res = (j.x - i.x) * (k.x - i.x) + (j.y - i.y) * (k.y - i.y);
        double len1 = sqrt((j.x - i.x) * (j.x - i.x) + (j.y - i.y) * (j.y - i.y));
        double len2 = sqrt((k.x - i.x) * (k.x - i.x) + (k.y - i.y) * (k.y - i.y));
        return fabs(res) / (len1 * len2);
    }
    void init(PointSet cur){
        seed = seed * seed % INT_MAX;
        srand(seed);
        city_num = cur.siz;
        for(int i = 0; i < city_num; ++i) c[i] = cur.p[i];
        for(int i = 0; i < city_num; ++i){
            for(int j = 0; j < city_num; ++j){
                Tau[i][j] = 1.5;
                dis[i][j] = dist(cur.p[i], cur.p[j]);
                Eta[i][j] = 10 / dis[i][j];
                for(int k = 0; k < city_num; ++k){
                    if(i == j || i == k) continue;
                    Relation[i][j][k] = cossim(tp[i], tp[j], tp[k]);
                }
            }
        }
    }
    int rand_chose(int n){
        double res = (double)rand() / RAND_MAX;
        for(int i = 0; i < n; ++i){
            res -= P[i];
            if(res <= 0) return i;
        }
        return 0;
    }
    void compute_paths(){
        for(int i = 0; i < ant_num; ++i){
            double sum = dis[0][city_num - 1];
            for(int j = 0; j < city_num - 1; ++j){
                int a = AntPath[i][j], b = AntPath[i][j + 1];
                sum += dis[a][b];
            }
            pathlen[i] = sum;
        }
    }
    void update_Tau(){
        for(int i = 0; i < ant_num; ++i){
            for(int j = 0; j < city_num; ++j){
                DeltaTau[i][j] = 0;
            }
        }
        for(int i = 0; i < ant_num; ++i){
            for(int j = 0; j < city_num - 1; ++j){
                int a = AntPath[i][j], b = AntPath[i][j + 1];
                DeltaTau[a][b] += Q / pathlen[i];
            }
            int a = AntPath[i][city_num - 1], b = AntPath[i][0];
            DeltaTau[a][b] += Q / pathlen[i];
        }
        for(int i = 0; i < city_num; ++i){
            for(int j = 0; j < city_num; ++j){
                for(int k = 0; k < city_num; ++k){
                    DeltaTau[i][k] += Tau[i][j] * Relation[i][j][k] * delta;
                }
            }
        }
        for(int i = 0; i < city_num; ++i){
            for(int j = 0; j < city_num; ++j){
                Tau[i][j] = (1 - rho) * Tau[i][j] + DeltaTau[i][j];
            }
        }
    }
    void GetAnts(){
        for(int k = 0; k < ant_num; ++k){
            int s = rand() % city_num;
            AntPath[k][0] = s;
            set<int>unvisit;
            for(int i = 0; i < city_num; ++i){
                if(i == s) continue;
                unvisit.insert(i);
            }
            int pos = s, cnt = 0;
            AntPath[k][0] = pos;
            while(unvisit.size()){
                int i = 0;
                double sumP = 0;
                for(auto v:unvisit){
                    P[i] = pow(Tau[pos][v], alpha) * pow(Eta[pos][v], beta);
                    sumP += P[i];
                    i++;
                }
                for(int j = 0; j < i; ++j) P[j] /= sumP;
                int index = rand_chose(i);
                int nxt;
                i = 0;
                for(auto v:unvisit){
                    if(i == index) nxt = v;
                    i++;
                }
                AntPath[k][++cnt] = nxt;
                unvisit.erase(nxt);
            }
        }
    }
    vector<int> run(PointSet cur){
        init(cur);  
        double best_length = INT_MAX;
        vector<int>best_path(city_num);
        for(int iter = 0; iter < iter_max; ++iter){
            GetAnts();
            compute_paths();
            int pos = 0;
            for(int i = 1; i < ant_num; ++i){
                if(pathlen[i] < pathlen[pos]) pos = i;
            }
            if(pathlen[pos] < best_length){
                best_length = pathlen[pos];
                for(int i = 0; i < city_num; ++i) best_path[i] = c[AntPath[pos][i]];
            }
            update_Tau();
            // cout << "iter: " << iter << ", len: " << pathlen[pos] << "\n";
        }
        return best_path;
    }
}

int setnum;
double minlen = INT_MAX;
bool vis[N];

void Dfs(int x, int num, double len){
    if(num > setnum){
        if(len < minlen) minlen = len;
        return;
    }
    for(int i = 0; i < setnum; ++i){
        if(vis[i]) continue;
        vis[i] = true;
        Dfs(i, num + 1, len + ps[x].sum + dist(ps[x].brk1, ps[i].brk1));
        Dfs(i, num + 1, len + ps[x].sum + dist(ps[x].brk2, ps[i].brk1));
        Dfs(i, num + 1, len + ps[x].sum + dist(ps[i].brk2, ps[x].brk1));
        Dfs(i, num + 1, len + ps[x].sum + dist(ps[x].brk2, ps[i].brk2));
    }
}
int main(){
    freopen("dantzig42.tsp", "r", stdin);
    // freopen("ans.txt", "w", stdout);
    int n;
    cin >> n;
    for(int i = 1; i <= n; ++i){
        int id;
        double x, y;
        cin >> id >> x >> y;
        tp[i].x = x, tp[i].y = y;
        lst.insert(i);
    }
    int s = sqrt(n);
    while(lst.size()){
        auto it = lst.begin();
        int a = *it;
        lst.erase(it);
        vector<pair<double, int>>dis;
        for(auto b:lst){
            double cur = dist(a, b);
            dis.push_back({cur, b});
        }
        sort(dis.begin(), dis.end(), [&](pair<double, int> x, pair<double, int> y){
            return x.first < y.first;
        });
        PointSet cur;
        cur.siz = 0;
        ++cur.siz;
        cur.p.push_back(a);
        // cout<<a<<" ";
        for(int i = 0; i < s && i < dis.size(); ++i){
            ++cur.siz;
            int x = dis[i].second;
            cur.p.push_back(x);
            lst.erase(x);
            // cout<<x<<" ";
        }
        // cout<<"\n";
        ps.push_back(cur);
    }
    setnum = ps.size();
    int T = 20;
    double maxans = 0, minans = INT_MAX, ave = 0;
    while(T--){
        for(int i = 0; i < setnum; ++i){
            ps[i].sum = 0;
            vector<int>pth;
            pth = ACO::run(ps[i]);
            ps[i].subpth = pth;
            double maxx = 0;
            for(int j = 0; j < ps[i].siz; ++j){
                double res = dist(pth[j], pth[j + 1]);
                if(i > 0 || j > 0) ps[i].sum += res;
                if(res > maxx) maxx = res, ps[i].brk1 = pth[j], ps[i].brk2 = pth[j + 1];
            }
            ps[i].sum -= maxx;
        }
        minlen = INT_MAX;
        for(int i = 0; i < setnum; ++i) vis[i] = false;
        Dfs(0, 1, 0);
        ave += minlen;
        if(minlen > maxans) maxans = minlen;
        if(minlen < minans) minans = minlen;
    }
    ave /= 20;
    cout << "ave_dis = " << ave << " max_dis = " << maxans << " mindis = " << minans << "\n";
    return 0;
}