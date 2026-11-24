// Input: the nine lists of weight-4 and weight-8 words in codes.txt.
// Traverse each six-set first, then its twenty three-subsets.
// Output every cut surviving the displayed diagonal trace inequalities.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <set>
#include <string>
#include <vector>
using Word=uint32_t;
int wt(Word x) { return __builtin_popcount(x); }
const Word end=1U<<24;
Word next_popcount(Word x) {
    Word low=x&-x, raised=x+low;
    return raised+((raised^x)/low>>2);
}
static bool admits(const std::vector<Word>& supports) {
    if (supports.size()<12) return false;
    Word meet=end-1, join=0;
    std::array<int,24> incidence{};
    for (auto s:supports) {
        if(wt(s)==1) return false;
        assert(wt(s)==5);
        meet &= s; join |= s;
        for(Word u=s;u;u&=u-1) ++incidence[__builtin_ctz(u)];
    }
    if(meet || wt(join)!=18) return false;
    std::array<int,24> columns{};
    std::iota(columns.begin(),columns.end(),0);
    std::stable_sort(columns.begin(),columns.end(),[&](int x,int y){return incidence[x]>incidence[y];});
    Word subset=0;
    std::vector<int> levels(supports.size());
    for(int r=1;r<18;++r) {
        subset|=1U<<columns[r-1];
        std::transform(supports.begin(),supports.end(),levels.begin(),[&](Word s){return wt(s&subset);});
        std::sort(levels.begin(),levels.end());
        long long lo=0,hi=0;
        for(int j=0;j<58;++j) {
            lo+=levels[j/5]; hi+=levels[levels.size()-1-j/5];
        }
        long long gap=std::max({0LL,18*lo-290*r,290*r-18*hi});
        if(gap*gap>1508LL*r*(18-r))return false;
    }
    return true;
}
int main() {
    std::string name; size_t n;
    while(std::cin>>name>>n) {
        std::vector<Word> words(n);
        for(Word& w:words)std::cin>>w;
        std::set<Word> whole{0};
        for(auto w:words) if(!whole.count(w)) {
            auto snapshot=whole;
            for(auto x:snapshot)whole.insert(w^x);
        }
        if(whole.size()!=4096) {
            // The odd half of the paired d24 code has constant weight twelve.
            assert(name=="d24" && whole.size()==2048);
            const Word transversal=0x555555;
            auto snapshot=whole;
            for(auto x:snapshot) {
                assert(wt(x^transversal)==12);
                whole.insert(x^transversal);
            }
        }
        assert(whole.size()==4096);
        std::vector<Word> colored;
        for(auto w:whole) {
            assert(wt(w)%4==0);
            if(wt(w)==4 || wt(w)==8)colored.push_back(w);
        }
        auto copy=words;std::sort(copy.begin(),copy.end()); assert(copy==colored);
        std::vector<int> slot(end,-1);
        unsigned examined=0,kept=0;
        for(Word six=63;six<end;six=next_popcount(six)) {
            std::array<std::vector<Word>,20> blocks;
            std::array<Word,20> triples{};int m=0;
            for(Word sub=six;sub;sub=(sub-1)&six)if(wt(sub)==3){slot[sub]=m;triples[m++]=sub;}
            assert(m==20);
            for(auto w:words){Word cut=w&six;if(wt(cut)==3)blocks[slot[cut]].push_back(w^cut);}
            for(int i=0;i<m;++i)if(admits(blocks[i])){
                std::cout<<name<<" "<<triples[i]<<" "<<(six^triples[i])<<"\n";++kept;
            }
            examined+=m;
        }
        assert(examined==2691920);
        std::cerr<<name<<" "<<examined<<" "<<kept<<"\n";
    }
}
