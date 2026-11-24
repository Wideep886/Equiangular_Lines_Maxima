#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>
using Row=std::array<int,24>;
using Mask=std::uint32_t;
void require(bool ok){if(!ok)throw std::runtime_error("inconsistent shell enumeration");}
int main(){
 int nr;
 while(std::cin>>nr){
  Row h;for(int& x:h)std::cin>>x;
  std::array<std::vector<int>,24> nb;
  int edges;std::cin>>edges;
  for(int k=0;k<edges;k++){int i,j;std::cin>>i>>j;nb.at(i).push_back(j);nb.at(j).push_back(i);}
  std::map<Row,int> index;
  std::vector<Row> rows(nr),adj(nr);
  for(int k=0;k<nr;k++){for(int& x:rows[k])std::cin>>x;require(index.emplace(rows[k],k).second);}
  Mask rootSet=0;
  for(int i=0;i<24;i++)if(!h[i])rootSet|=Mask(1)<<i;
  std::vector<Mask> neg(nr),nz(nr);
  for(int k=0;k<nr;k++)for(int i=0;i<24;i++){
   int label=2*rows[k][i];for(int j:nb[i])label-=rows[k][j];
   if(label<0)neg[k]|=Mask(1)<<i;
   if(label!=0)nz[k]|=Mask(1)<<i;
   if(!(rootSet & (Mask(1)<<i))) continue;
   Row t=rows[k];t[i]-=label;adj[k][i]=index.at(t);
  }
  auto bound=[&](Mask active,Mask mark){
   std::vector<int> distance(nr,-1), queue,layer,gens;
   for(int i=0;i<24;i++)if((active>>i)&1)gens.push_back(i);
   int sumBound=0;
   for(int k=0;k<nr;k++){
    if((nz[k]&mark)||(neg[k]&active))continue;
    require(distance[k]==-1); // unique dominant element per orbit
    distance[k]=0;queue.assign(1,k);layer.assign(1,1);
    for(std::size_t p=0;p<queue.size();p++){
     int j=queue[p];
     for(int i:gens){
      int t=adj[j][i];require(!(nz[t]&mark));
      if(distance[t]!=-1)continue;
      distance[t]=1+distance[j];
      if(int(layer.size())==distance[t])layer.push_back(0);
      layer[distance[t]]++;queue.push_back(t);
     }
    }
    int mid=layer[layer.size()/2];
    require(mid==*std::max_element(layer.begin(),layer.end()));
    require(std::equal(layer.begin(),layer.end(),layer.rbegin()));
    sumBound+=mid;
   }
   for(int k=0;k<nr;k++)require(bool(nz[k]&mark)==(distance[k]==-1));
   return sumBound;
  };
  int first=bound(rootSet,0);std::cout<<first;
  if(first>=116)for(int i=0;i<24;i++)for(int j:nb[i])if(i>j){
   Mask m=(Mask(1)<<i)|(Mask(1)<<j);
   if((rootSet&m)!=m)continue;
   Mask forbidden=m;for(int a:nb[i])forbidden|=Mask(1)<<a;for(int a:nb[j])forbidden|=Mask(1)<<a;
   std::cout<<' '<<i<<','<<j<<','<<bound(rootSet&~forbidden,m);
  }
  std::cout<<'\n';
 }
}
