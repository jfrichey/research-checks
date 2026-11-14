// Site-percolation frontier chain. A dyadic approximate stationary measure is
// certified by its EXACT residual and the all-closed-row Doeblin minorization.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;
using i128=__int128_t;
string decimal(i128 a) {
    if(a==0) return "0";
    bool negative=a<0;
    if(negative) a=-a;
    string s;
    while(a){s+=char('0'+a%10);a/=10;}
    if(negative)s+='-';
    reverse(s.begin(),s.end());return s;
}
struct DSU {
    int p[32];
    DSU(int n){for(int i=0;i<n;i++)p[i]=i;}
    int root(int i){while(i!=p[i]){p[i]=p[p[i]];i=p[i];}return i;}
    int join(int i,int j){i=root(i);j=root(j);if(i==j)return 0;p[i]=j;return 1;}
};
pair<uint64_t,int> step(uint64_t state,unsigned mask,int w) {
    DSU d(2*w);
    int seen[16],a[16];fill(seen,seen+16,-1);
    for(int i=0;i<w;i++) {
        a[i]=(state>>(4*i))&15;
        if(a[i]){if(seen[a[i]]>=0)d.join(i,seen[a[i]]);else seen[a[i]]=i;}
    }
    int reward=__builtin_popcount(mask);
    for(int i=0;i<w;i++)if((mask>>i)&1) {
        if(a[i])reward-=d.join(i,w+i);
        if(i && ((mask>>(i-1))&1))reward-=d.join(w+i-1,w+i);
    }
    int labels[32]={},next=0;uint64_t out=0;
    for(int i=0;i<w;i++)if((mask>>i)&1) {
        int r=d.root(w+i);if(!labels[r])labels[r]=++next;
        out|=uint64_t(labels[r])<<(4*i);
    }
    return {out,reward};
}
int main(int argc,char**argv) {
    if(argc<3){cerr<<"Usage: certify_strips WIDTH WIRED [CERTIFICATE_PATH]\n";return 2;}
    int k=stoi(argv[1]);bool wired=stoi(argv[2]);int w=k+wired,D=1<<k;
    if(k<1||k>12)return 2;
    vector<uint64_t> states{uint64_t(wired)};
    unordered_map<uint64_t,int> index{{states[0],0}};
    vector<uint32_t> transition;
    vector<int> rewards;
    for(size_t i=0;i<states.size();i++) {
        uint64_t old=states[i];int reward=0;
        for(int mask=0;mask<D;mask++) {
            auto [s,r]=step(old,wired?(2*mask+1):mask,w);
            auto [it,inserted]=index.emplace(s,states.size());
            if(inserted)states.push_back(s);
            transition.push_back(it->second);reward+=r;
        }
        rewards.push_back(reward);
    }
    const size_t n=states.size();
    cerr<<"width="<<k<<" wired="<<wired<<" states="<<n<<"\n";
    vector<long double> pi(n),next(n);pi[0]=1;
    int iterations=0;
    for(;iterations<1000;iterations++) {
        fill(next.begin(),next.end(),0);
        for(size_t i=0;i<n;i++) {
            long double a=pi[i]/D;
            for(int m=0;m<D;m++)next[transition[i*D+m]]+=a;
        }
        long double diff=0,mass=0;
        for(size_t i=0;i<n;i++){mass+=next[i];diff+=abs(next[i]-pi[i]);}
        for(auto &a:next)a/=mass;
        pi.swap(next);
        if(diff<1e-18L)break;
    }
    const int64_t Q=int64_t(1)<<60;
    vector<int64_t> q(n);int64_t total=0;
    for(size_t i=0;i<n;i++){q[i]=int64_t(pi[i]*Q);total+=q[i];}
    q[0]+=Q-total;
    if(q[0]<0){cerr<<"negative mass\n";return 3;}
    vector<i128> residual(n);
    i128 expectation=0;
    for(size_t i=0;i<n;i++) {
        expectation+=i128(q[i])*rewards[i];
        residual[i]-=i128(q[i])*D;
        for(int m=0;m<D;m++)residual[transition[i*D+m]]+=q[i];
    }
    i128 norm=0;for(auto r:residual)norm+=(r<0?-r:r);
    // At least probability 1/D of the SAME reset state in every row.
    for(size_t i=0;i<n;i++)if(transition[i*D]!=0)return 4;
    // |reward(state)| <= w, so error <= w * residual_numerator / Q.
    for(int r:rewards)if(abs(r)>w*D)return 5;
    if(argc>3){ofstream out(argv[3]);for(size_t i=0;i<n;i++)out<<states[i]<<" "<<q[i]<<"\n";}
    cout<<"{\"width\":"<<k<<",\"wired\":"<<(wired?"true":"false")
        <<",\"states\":"<<n<<",\"iterations\":"<<iterations+1
        <<",\"Q\":\""<<Q<<"\",\"D\":"<<D
        <<",\"expectation_numerator\":\""<<decimal(expectation)
        <<"\",\"residual_numerator\":\""<<decimal(norm)
        <<"\",\"reward_bound\":"<<w<<"}\n";
}
