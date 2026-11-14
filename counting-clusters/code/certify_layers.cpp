// Single-site square-strip transfer, with exact integer residual certificates.
// Noncrossing frontiers use three bits per position. Layer maps are temporary.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;
using I=__int128_t;
string dec(I x){
    if(!x)return "0";
    bool neg=x<0;if(neg)x=-x;string s;
    while(x){s+=char('0'+x%10);x/=10;}
    if(neg)s+='-';
    reverse(s.begin(),s.end());return s;
}
I parse(const string&s){
    if(s.empty()||s[0]=='-')throw runtime_error("invalid nonnegative mass");
    I x=0;for(char c:s){
        if(c<'0'||c>'9'||x>(I(1)<<100))throw runtime_error("invalid mass");
        x=10*x+c-'0';
    }return x;
}
// Symbols: 0 vacant; 1 singleton; 2 first; 3 middle; 4 last in a block.
// The permanent wire is an implicit, always-open outermost block, label 1.
void decode(uint64_t s,int*a,int w,bool wired){
    int stack[24],depth=0,next=0;
    if(wired)stack[depth++]=++next;
    for(int j=0;j<w;j++){
        int t=(s>>(3*j))&7;a[j]=0;
        if(t==1)a[j]=++next;
        else if(t==2){a[j]=++next;stack[depth++]=a[j];}
        else if(t==3||t==4){
            if(!depth)throw runtime_error("invalid frontier code");
            a[j]=stack[depth-1];if(t==4)--depth;
        }else if(t!=0)throw runtime_error("invalid frontier symbol");
    }
    if(depth!=int(wired))throw runtime_error("unclosed frontier code");
}
uint64_t encode(const int*a,int w,bool wired){
    int count[32]={},first[32],last[32],stack[24],depth=0;
    fill(first,first+32,-1);
    if(wired)stack[depth++]=1;
    for(int j=0;j<w;j++)if(a[j]){
        ++count[a[j]];if(first[a[j]]<0)first[a[j]]=j;last[a[j]]=j;
    }
    uint64_t s=0;
    for(int j=0;j<w;j++)if(a[j]){
        int x=a[j],t;
        if(wired&&x==1)t=3;
        else if(count[x]==1)t=1;
        else if(j==first[x])t=2;
        else if(j==last[x])t=4;
        else t=3;
        if(t==2)stack[depth++]=x;
        if(t==3||t==4){
            if(!depth||stack[depth-1]!=x)throw runtime_error("crossing frontier");
            if(t==4)--depth;
        }
        s|=uint64_t(t)<<(3*j);
    }
    if(depth!=int(wired))throw runtime_error("unclosed partition");
    return s;
}
pair<uint64_t,int> site(uint64_t s,int pos,bool occupied,int w,bool wired){
    int a[24]={};decode(s,a,w,wired);
    int above=a[pos],left=pos?a[pos-1]:int(wired),r=0;
    if(!occupied)a[pos]=0;
    else{
        r=1;a[pos]=above?above:31;if(above)--r;
        if(left&&left!=a[pos]){
            int from=a[pos],to=left;
            if(wired&&from==1)swap(from,to);
            for(int j=0;j<w;j++)if(a[j]==from)a[j]=to;
            --r;
        }
    }
    return {encode(a,w,wired),r};
}
vector<uint64_t> frontiers(int w,int kink,bool wired){
    vector<uint64_t> result;
    // Occupied neighbors agree except across the one step in the frontier.
    // This is an explicit superset of reachable planar states. Closure and
    // the common reset are checked, so no reachability conjecture is needed.
    auto visit=[&](auto&&self,int pos,int depth,int prev,uint64_t s)->void{
        if(depth-int(wired)>w-pos)return;
        if(pos==w){if(depth==int(wired))result.push_back(s);return;}
        self(self,pos+1,depth,0,s);
        bool forced=(pos==0)?wired:(pos!=kink&&prev!=0);
        if(!forced){
            self(self,pos+1,depth,1,s|(uint64_t(1)<<(3*pos)));
            self(self,pos+1,depth+1,2,s|(uint64_t(2)<<(3*pos)));
        }
        bool can_resume=depth>0&&(!forced||pos==0||prev==2||prev==3);
        if(can_resume){
            self(self,pos+1,depth,3,s|(uint64_t(3)<<(3*pos)));
            if(depth>int(wired))self(self,pos+1,depth-1,4,s|(uint64_t(4)<<(3*pos)));
        }
    };
    visit(visit,0,int(wired),0,0);return result;
}
// Independent whole-row union-find update, used only by small-width checks.
pair<uint64_t,int> direct_row(uint64_t s,unsigned mask,int k,bool wired){
    int w=k+int(wired),a[24]={},random[24]={};decode(s,random,k,wired);
    if(wired)a[0]=1;
    for(int j=0;j<k;j++)a[j+int(wired)]=random[j];
    if(wired)mask=2*mask+1;
    int p[48],seen[32];for(int i=0;i<2*w;i++)p[i]=i;
    fill(seen,seen+32,-1);
    auto root=[&](int i){while(p[i]!=i)i=p[i];return i;};
    auto join=[&](int i,int j){i=root(i);j=root(j);if(i==j)return 0;p[i]=j;return 1;};
    for(int i=0;i<w;i++)if(a[i]){
        if(seen[a[i]]>=0)join(i,seen[a[i]]);else seen[a[i]]=i;
    }
    int r=__builtin_popcount(mask);
    for(int i=0;i<w;i++)if(mask>>i&1){
        if(a[i])r-=join(i,w+i);
        if(i&&(mask>>(i-1)&1))r-=join(w+i-1,w+i);
    }
    int lab[48]={},next=0,out[24]={};
    for(int i=0;i<w;i++)if(mask>>i&1){
        int j=root(w+i);if(!lab[j])lab[j]=++next;
        if(i>=int(wired))out[i-int(wired)]=lab[j];
    }
    return {encode(out,k,wired),r};
}
struct Layer{vector<uint32_t> vacant,occupied;vector<int8_t> reward;};
int main(int argc,char**argv)try{
    if(argc<3)throw runtime_error("Usage: certify_layers WIDTH WIRED [--output PATH] [--verify PATH] [--check-rows]");
    int k=stoi(argv[1]);bool wired=stoi(argv[2]);
    if(k<1||k>17)throw runtime_error("width must be 1..17");
    int D=1<<k;string output,verify;bool check=false;
    for(int i=3;i<argc;i++){
        string arg=argv[i];
        if(arg=="--check-rows")check=true;
        else if(arg=="--output"&&i+1<argc)output=argv[++i];
        else if(arg=="--verify"&&i+1<argc)verify=argv[++i];
        else throw runtime_error("invalid argument");
    }
    vector<uint64_t> states=frontiers(k,0,wired),current=states;
    if(states[0]!=0)throw runtime_error("reset not first");
    size_t n=states.size(),layered=0;vector<Layer> layers(k);
    for(int pos=0;pos<k;pos++){
        vector<uint64_t> next=pos+1==k?states:frontiers(k,pos+1,wired);
        unordered_map<uint64_t,uint32_t> index;index.reserve(next.size());
        for(uint32_t j=0;j<next.size();j++)if(!index.emplace(next[j],j).second)throw runtime_error("duplicate state");
        auto&l=layers[pos];l.vacant.resize(current.size());l.occupied.resize(current.size());l.reward.resize(current.size());
        for(size_t j=0;j<current.size();j++){
            auto [s0,r0]=site(current[j],pos,false,k,wired);
            auto [s1,r1]=site(current[j],pos,true,k,wired);
            if(r0!=0||r1 < -1||r1>1)throw runtime_error("invalid site reward");
            l.vacant[j]=index.at(s0);l.occupied[j]=index.at(s1);l.reward[j]=r1;
        }
        layered+=current.size();current.swap(next);
    }
    vector<uint64_t>().swap(current);
    cerr<<"width="<<k<<" wired="<<wired<<" row states="<<n<<" layered states="<<layered<<"\n";
    uint64_t checked=0;
    for(uint32_t j=0;j<n;j++){
        uint32_t t=j;for(int pos=0;pos<k;pos++)t=layers[pos].vacant[t];
        if(t!=0)throw runtime_error("all-white row fails to reset");
        if(check){
            if(k>8)throw runtime_error("row checks restricted to width <=8");
            for(int mask=0;mask<D;mask++){
                t=j;int r=0;
                for(int pos=0;pos<k;pos++){
                    auto&l=layers[pos];bool b=mask>>pos&1;
                    if(b)r+=l.reward[t];
                    t=b?l.occupied[t]:l.vacant[t];
                }
                if(direct_row(states[j],mask,k,wired)!=pair<uint64_t,int>{states[t],r})throw runtime_error("whole-row cross-check failed");
                ++checked;
            }
        }
    }
    const I Q=I(1)<<80;vector<I> q(n);int iterations=0;
    if(verify.empty()){
        vector<long double> pi(n),v,next;pi[0]=1;
        for(;iterations<1500;iterations++){
            v=pi;
            for(int pos=0;pos<k;pos++){
                auto&l=layers[pos];next.assign(layers[(pos+1)%k].vacant.size(),0);
                for(size_t j=0;j<v.size();j++){
                    auto mass=v[j]/2;next[l.vacant[j]]+=mass;next[l.occupied[j]]+=mass;
                }v.swap(next);
            }
            long double mass=0,diff=0;for(auto x:v)mass+=x;
            for(size_t j=0;j<n;j++){v[j]/=mass;diff+=abs(v[j]-pi[j]);}
            pi.swap(v);if(diff<1e-18L)break;
        }
        I total=0;for(size_t j=0;j<n;j++){q[j]=I(pi[j]*static_cast<long double>(Q));total+=q[j];}
        q[0]+=Q-total;
    }else{
        unordered_map<uint64_t,uint32_t> index;index.reserve(n);
        for(uint32_t j=0;j<n;j++)index.emplace(states[j],j);
        ifstream in(verify);if(!in)throw runtime_error("certificate unavailable");
        vector<bool> seen(n);uint64_t state;string mass;size_t count=0;
        while(in>>state>>mass){
            auto it=index.find(state);
            if(it==index.end()||seen[it->second])throw runtime_error("unknown or duplicate state");
            seen[it->second]=true;q[it->second]=parse(mass);++count;
        }
        if(count!=n)throw runtime_error("missing certificate states");
    }
    I total=0;for(auto mass:q){if(mass<0||mass>Q)throw runtime_error("invalid candidate mass");total+=mass;}
    if(total!=Q)throw runtime_error("incorrect candidate normalization");
    vector<I> v=q,next;I expectation=0;
    for(int pos=0;pos<k;pos++){
        auto&l=layers[pos];next.assign(layers[(pos+1)%k].vacant.size(),0);expectation*=2;
        for(size_t j=0;j<v.size();j++){
            next[l.vacant[j]]+=v[j];next[l.occupied[j]]+=v[j];expectation+=v[j]*l.reward[j];
        }v.swap(next);
    }
    I residual=0,mass=0;
    for(size_t j=0;j<n;j++){mass+=v[j];I r=v[j]-D*q[j];residual+=r<0?-r:r;}
    if(mass!=Q*D)throw runtime_error("integer transfer mass mismatch");
    if(!output.empty()){
        ofstream out(output);if(!out)throw runtime_error("cannot write certificate");
        for(size_t j=0;j<n;j++)out<<states[j]<<' '<<dec(q[j])<<'\n';
    }
    cout<<"{\"width\":"<<k<<",\"wired\":"<<(wired?"true":"false")<<",\"states\":"<<n
        <<",\"layered_states\":"<<layered<<",\"iterations\":"<<iterations<<",\"Q\":\""<<dec(Q)<<"\",\"D\":"<<D
        <<",\"expectation_numerator\":\""<<dec(expectation)<<"\",\"residual_numerator\":\""<<dec(residual)
        <<"\",\"reward_bound\":"<<k<<",\"whole_row_checks\":"<<checked<<",\"loaded_certificate\":"<<(verify.empty()?"false":"true")
        <<",\"state_encoding\":\"noncrossing-3bit-implicit-wire\"}\n";
    return 0;
}catch(const exception&e){cerr<<e.what()<<'\n';return 1;}
