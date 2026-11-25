// Independently traverse all local second-label tuples and all admissible signed residual pairings.
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <set>
#include <vector>
using Row=std::array<int,12>;
using Local=std::vector<std::vector<Row>>;
struct Search {
    int S,A,n;
    const std::vector<Local>& components;
    std::vector<int> anchor;
    std::set<int> found;
    std::vector<int> squarefree;
    static constexpr std::array<std::array<int,2>,10> coefficients{{
        {0,2},{0,3},{1,1},{1,-1},{1,2},{1,-2},{2,1},{2,-1},{3,1},{3,-1}
    }};
    bool valid(int z, int parity) const {
        return (z==0 || z>=5*S) && (z%(4*S))==parity*S;
    }
    void inspect(const Row& row) {
        const int B=5*S-row[0];
        if (B!=0 && row[1]>3) return;
        if(!valid(row[2]+4*B,0) || !valid(row[3]+9*B,1))return;
        int radius=0;
        while((radius+1)*(radius+1)<=4*A*B)++radius;
        // Determine every signed integer y with positive area, directly from Cauchy.
        // Modulo 4S from the (1,1) pair only accelerates the loop.
        int lo=-radius;
        while ((row[4]+A+B+lo)%(4*S)!=0)++lo;
        for(int y=lo;y<=radius;y+=4*S) {
            const int K=4*A*B-y*y;
            if(K<=0)continue;
            bool ok=true;
            for(int j=2;j<10;j++) {
                auto [m,k]=coefficients[j];
                const int length=row[j+2]+m*m*A+k*k*B+m*k*y;
                if(!valid(length,(std::abs(m+k)%2))){ok=false;break;}
            }
            if(ok)found.insert(squarefree[K]);
        }
    }
    void walk(int i,Row row){
        if(i==n){inspect(row);return;}
        for(const Row& v:components[i][anchor[i]]) {
            if(row[0]+v[0]>5*S)continue;
            Row next=row;
            for(int j=0;j<12;j++)next[j]+=v[j];
            walk(i+1,next);
        }
    }
};
int main(){
    std::ios::sync_with_stdio(false);std::cin.tie(nullptr);
    unsigned long long anchors=0,total=0;
    int ntypes;std::cin>>ntypes;
    for(int j=0;j<ntypes;j++) {
        int n,S;std::cin>>n>>S;assert(S<=210);
        std::vector<Local> parts(n);
        for(auto& local:parts){
            int nlab;std::cin>>nlab;local.resize(nlab);
            for(auto& rows:local){
                int d;std::cin>>d;rows.resize(d);
                for(auto& row:rows)for(auto& a:row)std::cin>>a;
            }
        }
        Search search{S,0,n,parts,std::vector<int>(n),{},{}};
        const int maxK=100*S*S;
        search.squarefree.resize(maxK+1);
        for(int v=0;v<=maxK;v++)search.squarefree[v]=v;
        // Factoring out all square divisors by increasing integer squares.
        for(int t=2;t*t<=maxK;t++)
            for(int v=t*t;v<=maxK;v+=t*t)search.squarefree[v]=v/(t*t);
        int count;std::cin>>count;
        while(count--) {
            std::cin>>search.A;for(int& k:search.anchor)std::cin>>k;
            int size;std::cin>>size;std::set<int> expect;
            for(int h=0;h<size;h++){int k;std::cin>>k;expect.insert(k);}
            search.found.clear(); search.walk(0,Row{});
            if(search.found!=expect){
                std::cerr<<"different sectors, type "<<j<<" anchor "<<anchors<<" actual "<<search.found.size()<<" expected "<<expect.size()<<"\n";return 1;
            }
            anchors++;total+=size;
        }
    }
    assert(std::cin);
    std::cout<<ntypes<<" types "<<anchors<<" positive residual anchors "<<total<<" labelled area sectors verified\n";
}
