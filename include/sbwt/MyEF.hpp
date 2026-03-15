
#pragma once

#include <algorithm>
#include <cstdio>
#include <cstdint>
#include <iostream>
#include <chrono>

#include <algorithm>
#include <iostream>
#include <random>


using namespace std::chrono;

//long long cur_time_millis(){
//	return (std::chrono::duration_cast< milliseconds >(high_resolution_clock::now().time_since_epoch())).count();
//}
//
//long long cur_time_micros(){
//	return (std::chrono::duration_cast< microseconds >(high_resolution_clock::now().time_since_epoch())).count();
//}

#define MPER 8

class MyEF{
   uint8_t *_L; //lower bits
   uint64_t *_U; //buckets sizes, in unary
   uint32_t _n; //number of elements
   uint64_t _u; //universe
   uint64_t *_M; //support for select on _U
   uint64_t _mlen; //the length of _M
   uint64_t _mper; //a period used to build _M FIXME: hard code to 8

   uint64_t _nbuckets; //ceil(_n/16)

   uint64_t *_safe_select_U; //for debugging
   uint64_t *_safe_rank; //for debugging

   void setbit(uint64_t *bv, uint64_t i){
      //std::cerr << "i: " << i << '\n';
      //std::cerr << "imod64: " << (i%(uint64_t)64) << '\n';
      //fprintf(stderr,"Before: %llx\n",bv[i/64]);
      //std::cerr << bv[i/(uint64_t)64] << '\n';
      //std::cerr << ((uint64_t)1 << (i%(uint64_t)64)) << '\n';
      bv[i/(uint64_t)64] |= ((uint64_t)1 << (i%(uint64_t)64));
      //fprintf(stderr,"After : %llx\n",bv[i/64]);
   }

   uint64_t getbit(uint64_t *bv, uint64_t i) const{
      return ((bv[i/(uint64_t)64] >> (i%(uint64_t)64)) & (uint64_t)1);
   }

   inline void assignLBits(uint64_t index, uint64_t bits){
      //std::cerr << "index: "<<index<<", bits: "<<bits<<'\n';
      _L[index/2] = _L[index/2] | (bits << (4*(index%2)));
   }
   
   inline uint64_t getLBits(uint64_t index) const{
      return (_L[index/2]>>(4*(index%2)))&15;
   }
   
   void encodeCountsUnary(uint64_t *bucketCounts, uint64_t nbuckets, uint64_t n){
      uint64_t nEmptyBuckets = 0;
      for(uint64_t i=0; i<_nbuckets; i++){
         if(!bucketCounts[i]) nEmptyBuckets++;
      }
      std::cerr << "nEmptyBuckets: "<<nEmptyBuckets<<'\n';
      std::cerr << "n: " << n << '\n';
      std::cerr << "_nbuckets: " << _nbuckets << '\n';

      _safe_select_U = new uint64_t[_nbuckets];

      uint64_t ulen = n + _nbuckets; 
      uint64_t ulenwords = (ulen/64) + ((ulen%64)>0);
      std::cerr << "ulen: " << ulen << '\n';
      std::cerr << "ulenwords: " << ulenwords << '\n';
      _U = new uint64_t[ulenwords];
      for(uint64_t i=0;i<ulenwords;i++){
         _U[i] = 0;
      }
      uint64_t bi = 0;
      uint64_t uwi = 0;
      uint64_t nUbitsBefore = 0;
      for(uint32_t i=0; i<_nbuckets; i++){
         for(uint32_t j=0;j<bucketCounts[i];j++){
            bi++;
         }
         setbit(_U,bi);
         _safe_select_U[i] = bi;
         bi++;
      }

      _mper = MPER; //every 8 words we will store how many bits come before that word in _U
      _mlen = ulenwords/MPER; //_mper;
      //std::cerr << "_mlen: " << _mlen << '\n';
      _M = new uint64_t[_mlen];
      for(uint64_t i = 0; i < _mlen; i++){
         _M[i] = 0;
      }
      uint64_t bitcount = 0;
      for(uint64_t i=0;i<ulenwords;i++){
         if((i%MPER) == 0){
            _M[i/MPER] = bitcount;
            //std::cerr << "Setting bitcount: " << bitcount << '\n';
         }
         bitcount += __builtin_popcountll(_U[i]);
         //std::cerr << "bitcount: " << bitcount << '\n';
      }
   }

   inline uint64_t getNext1Pos(uint64_t prev1pos){
      uint64_t wi = prev1pos/(uint64_t)64;
      uint64_t w = _U[wi] >> (prev1pos%(uint64_t)64); //lowest bit now corresponds to prev1pos
      uint64_t n1p = 0;
      if(w >> (uint64_t)1){
         //after removing the bit corresponding to prev1pos, there are 1 bits left in w
         n1p = __builtin_ctz(w >> 1) + 1 + prev1pos;
      }else{
         wi++;
         n1p = wi*64; //prev1pos + (64-(prev1pos%64));
         w = _U[wi];
         //Because of our very specific implementation, the next five lines are not needed
         //--- the word *must* contain a 1 bit
         //while(!w){
         //   n1p += 64;
         //   wi++;
         //   w = _U[wi];
         //}
         n1p += __builtin_ctz(w);
      }
      return n1p;
   }

   inline uint64_t getPrev1Pos(uint64_t curr1pos) const{
      uint64_t wi = curr1pos/(uint64_t)64;
      //std::cerr << "wi: " << wi << '\n';
      uint64_t w = _U[wi];
      //std::cerr << "w: " << w << '\n';
      w = w << (((uint64_t)63)-(curr1pos%(uint64_t)64)); //highest bit now corresponds to curr1pos
      //std::cerr << "w: " << w << '\n';
      uint64_t p1p = 0;
      if(w << (uint64_t)1){
         //after removing the bit corresponding to prev1pos, there are 1 bits left in w
         p1p = curr1pos - (__builtin_clzll(w << (uint64_t)1) + (uint64_t)1);
         //std::cerr << "1: curr1pos p1p: "<<curr1pos<<' '<<p1p<<'\n';
      }else{
         p1p = curr1pos - (__builtin_ctzll(_U[wi]) + __builtin_clzll(_U[wi-(uint64_t)1]) + (uint64_t)1);
         //std::cerr << "2: curr1pos p1p: "<<curr1pos<<' '<<p1p<<'\n';
      }
      return p1p;
   }

public:

   
   //A is n integers drawn from universe u
   //A is sorted in ascending order and contains no duplicates
   void build(uint64_t *A, uint64_t n, uint64_t u){
      std::cerr << "build() called\n";
      _safe_rank = new uint64_t[u];
      for(uint64_t i=0,j=0; i<u; i++){
         _safe_rank[i] = j;
         if(A[j] == i){
            j++;
         }
      }
      _n = n;
      _u = u;
      _nbuckets = (u/(uint64_t)16)+((u%(uint64_t)16)>0);
      uint64_t *bucketCounts = new uint64_t[_nbuckets];
      _L = new uint8_t[(n/2) + (n%2)];
      for(uint64_t i = 0; i < _nbuckets; i++){
         bucketCounts[i] = 0;
      }
      for(uint64_t i = 0; i < n; i++){
         bucketCounts[A[i]/(uint64_t)16]++;
         assignLBits(i,A[i]%(uint64_t)16);
      }
      
      encodeCountsUnary(bucketCounts,_nbuckets,n);
      //for(uint32_t i=0;i<(_nbuckets+n);i++){
      //   std::cerr << getbit(_U,i);
      //}
      //std::cerr << '\n';
      //rank(124);
      //rank(324);
      //rank(624);
      //rank(1624);
      //rank(11624);
      //rank(111624);
      //for(int i=0;i<100000;i++) rank(i);
      //for(uint32_t i=0; i<n; i++){
      //   std::cerr << "ZZZ "<< A[i] << '\n';
      //}
      //rank(8928688);
      delete [] bucketCounts;


      std::cerr << "build() finished\n";
   }

   uint64_t enumerate(uint64_t *A){
      uint32_t prev1pos = 0;
      uint32_t next1pos = 0;
      uint32_t ndecoded = 0;
      uint32_t inc = getbit(_U,0);
      uint32_t bucki = inc;
      uint64_t x = inc;
      for(uint32_t i = 0; i < _n;){
         //one iteration of this loop processes one bucket (possibly empty)
         //1. Find the end of the bucket, which is the next 1 in _U
         next1pos = getNext1Pos(prev1pos);
         prev1pos+=inc; //move over the 1 (unless this is the first group and there are items in the first group)
         //std::cerr << "safe_select: " << _safe_select_U[x++] << '\n'; 
         //std::cerr << "prev1pos: " << prev1pos << '\n';
         //std::cerr << "next1pos: " << next1pos << '\n';
         //getchar();
         while(prev1pos < next1pos){
            uint32_t ell = getLBits(i);
            //std::cerr << "ell: " << ell << '\n';
            uint32_t v = bucki*16 + getLBits(i);
            //std::cerr << "v: " << v << '\n';
            A[ndecoded++] = v;
            prev1pos++;
            i++;
         }
         bucki++;
         inc = 1;
         //prev1pos = next1pos;
      }
      std::cerr << "Decoded "<<ndecoded<<"values.\n";
      return ndecoded;
   }


   void rank_test(uint64_t nqueries){
       std::random_device rd; // obtain a random number from hardware
       std::mt19937 gen(rd()); // seed the generator
       std::uniform_int_distribution<> distr(0, _u); // define the range
   
       uint64_t *Q = new uint64_t[nqueries];
   
       for(uint64_t i=0; i<nqueries; i++){
          Q[i] = distr(gen);
       }
       //int64_t t1 = cur_time_micros();
       uint64_t r = 0;
       for(uint64_t i=0; i<nqueries; i++){
          r += rank(Q[i]);
       }
       //int64_t t2 = cur_time_micros();
       //std::cerr << "Rank test time: "<<(double)(t2-t1)/(double)nqueries<<" micros per query\n";
       std::cerr << "Checksum: "<<r<<'\n';
       delete [] Q;
   }

   uint64_t lrank(uint64_t bstart, uint64_t bsize, uint64_t lkey) const{
      //std::cerr << "bstart bsize lkey: "<<bstart<<' '<<bsize<<' '<<lkey<<'\n';
      uint64_t r = 0;
      for(uint64_t i=0;i<bsize;i++){
         uint64_t lval = getLBits(bstart+i);
         //std::cerr << "lval: "<<lval<<'\n';
         if(lval >= lkey) break;
         r++;
      }
      return r;
   }

   uint64_t rank(uint64_t pos) const{
      //std::cerr << "---------------pos: "<<pos<<'\n';
      //let t = pos/16, where 16 is the bucket size
      uint64_t t = pos/16; //the bucket we seek
      //std::cerr << "t: "<<t<<'\n';
      //for(int i=0;i<10;i++) std::cerr << _M[i] << ' ';
      //std::cerr << '\n';
      //the first thing we need to determine for rank is where in _U the t^th 1 is
      //this is the same as select(_U,t) 
      //_M[x] is the number of 1s in _U before the (MPER*x)^th word 
      //so we start computing our answer to select(_U,t) by first computing the largest x s.t. _M[x] < t
      uint64_t *xp = std::lower_bound(_M,_M+_mlen,t); //an index into _M
      //std::cerr << "\txp: "<<xp<<'\n';
      //xp = _M+(t/320);
      //if(*xp < t){
      //   while(xp < _M+_mlen && *xp < t) xp++;
      //}else{
      //   while(xp > _M && *xp > t) xp--;
      //}
      //std::cerr << "\txp: "<<xp<<'\n';
      //std::cerr << "t -> index: "<<t<<' '<<(xp-_M)<<'\n';
      uint64_t x = (uint64_t)(xp - _M);
      if(x) x--;
      //std::cerr << "\tx: "<<x<<'\n';
      uint64_t xb = _M[x]; //the number of 1s before x*_mper 
      //std::cerr << "\txb: "<<xb<<'\n';
      x = x*MPER; //_mper; //index into _U of the word that has xb 1s before it
      //std::cerr << "\tx: "<<x<<'\n';

      //we are within MPER words of the bit we seek
      //next we move along word-at-a-time until we find the word containing our bit
      while(1){
         uint64_t bitsInWord = __builtin_popcountll(_U[x]);
         if((xb + bitsInWord) > t) break;
         xb += bitsInWord;
         x++;
      }
      //std::cerr << "\txb: "<<xb<<'\n';
      //std::cerr << "\tx: "<<x<<'\n';

      //x now points to the word containing the t^th 1 bit of _U
      //it is the (t-xb)^th 1 bit of _U[x]
      //find its position inside that word
      uint64_t w = _U[x];

      //std::cerr << "\tpc: " << __builtin_popcountll(_U[x])<<'\n';
      //for(int i=0;i<64;i++){
      //   std::cerr << (63-i)%10;
      //}
      //std::cerr << '\n';
      //for(int i=0;i<64;i++){
      //   std::cerr << getbit(_U,63-i);
      //}
      //std::cerr << '\n';

      uint64_t z = t-xb+1;
      //std::cerr << "\tz: "<<z<<'\n';
      
      int64_t y = -1;
      for(uint64_t i=0;i<z;i++){
         uint64_t shift = __builtin_ctzll(w) + 1;
         y = y + shift;
         w = w >> shift;
         //std::cerr << "\ty shift z: "<<y<<' '<<shift<<' '<<z<<'\n';
      }

      //std::cerr << "select("<<t<<"): "<<(64*x+y)<<'\n';
      //std::cerr << "_safe_select_U["<<t<<"]: "<<_safe_select_U[t]<<'\n';

      //if((64*x+y)!=_safe_select_U[t]){
      //   std::cerr << "Wap wap: problem with select (pos: "<<pos<<")\n";
      //   exit(1);
      //}


      uint64_t bs = 64*x+y;
      //std::cerr << ">>>>>bs: "<<bs<<'\n';
      if(t != 0){
         bs = bs - getPrev1Pos(bs) - 1; //get the number of items in the bucket
         //bs = getPrev1Pos(y) - y; //get the number of items in the bucket
      }
      //std::cerr << "bs: "<<bs<<'\n';
      uint64_t rank1 = (64*x + y) - t - bs; //the number of elements in buckets < t
      //now enumerate over the appropriate part of _Lbits
      uint64_t rank2 = lrank(rank1,bs,pos%16);

      //std::cerr << "rank1: "<<rank1<<'\n';
      //std::cerr << "rank2: "<<rank2<<'\n';
      //std::cerr << "rank("<<pos<<"): "<<rank1+rank2<<'\n';
      //std::cerr << "_safe_rank["<<pos<<"]: "<<_safe_rank[pos]<<'\n';
      //if(rank1+rank2 != _safe_rank[pos]){
      //   std::cerr << "Wap wap: problem with rank (pos: "<<pos<<")\n";
      //   exit(1);
      //}
      return rank1 + rank2;
   }

    // Returns the number of bytes written
    int64_t serialize(std::ostream& os) const{
        int64_t written = 0;
        uint64_t ulen = _n + _nbuckets;
        uint64_t ulenwords = (ulen/64) + ((ulen%64)>0);
        os.write((char *)&_n, sizeof(uint64_t));
        os.write((char *)&_u, sizeof(uint64_t));
        os.write((char *)&_mlen, sizeof(uint64_t));
        os.write((char *)&_mper, sizeof(uint64_t));
        os.write((char *)&_nbuckets, sizeof(uint64_t));
        std::cerr << "_n _u _mlen _mper _nbuckets: "<<_n<<' '<<_u<<' '<<_mlen<<' '<<_mper<<' '<<_nbuckets<<'\n';
        uint64_t llen = (_n/2) + (_n%2);
        os.write((char *)_L, llen);
        for(uint64_t i=0;i<10;i++) std::cerr << (uint8_t)_L[i] << ' ';
        std::cerr << '\n';
        os.write((char *)_U, ulenwords*sizeof(uint64_t));
        for(uint64_t i=0;i<10;i++) std::cerr << _U[i] << ' ';
        std::cerr << '\n';
        os.write((char *)_M, _mlen*sizeof(uint64_t));
        for(uint64_t i=0;i<10;i++) std::cerr << _M[i] << ' ';
        std::cerr << '\n';
        written += (sizeof(uint64_t)*(8));
        std::cerr << "MyEF serialize...\n"; 
        return written;
    }

    void load(std::istream& is){
        is.read((char *)&_n, sizeof(uint64_t));
        is.read((char *)&_u, sizeof(uint64_t));
        is.read((char *)&_mlen, sizeof(uint64_t));
        is.read((char *)&_mper, sizeof(uint64_t));
        is.read((char *)&_nbuckets, sizeof(uint64_t));
        std::cerr << "_n _u _mlen _mper _nbuckets: "<<_n<<' '<<_u<<' '<<_mlen<<' '<<_mper<<' '<<_nbuckets<<'\n';

        uint64_t llen = (_n/2) + (_n%2);
        _L = new uint8_t[llen];
        is.read((char *)_L, llen);
        for(uint64_t i=0;i<10;i++) std::cerr << (uint8_t)_L[i] << ' ';
        std::cerr << '\n';

        uint64_t ulen = _n + _nbuckets;
        uint64_t ulenwords = (ulen/64) + ((ulen%64)>0);
        _U = new uint64_t[ulenwords];
        is.read((char *)_U, ulenwords*sizeof(uint64_t));
        for(uint64_t i=0;i<10;i++) std::cerr << _U[i] << ' ';
        std::cerr << '\n';

        _M = new uint64_t[_mlen];
        is.read((char *)_M, _mlen*sizeof(uint64_t));
        for(uint64_t i=0;i<10;i++) std::cerr << _M[i] << ' ';
        std::cerr << '\n';
    }

    MyEF(){
    }

    MyEF(uint64_t *A, uint64_t n, uint64_t u){
       build(A,n,u);
    }

    MyEF(MyEF &other){
       this->_L = other._L;
       this->_U = other._U;
       this->_n = other._n;
       this->_u = other._u;
       this->_M = other._M;
       this->_mlen = other._mlen;
       this->_mper = other._mper;
       this->_nbuckets = other._nbuckets;
       this->_safe_select_U = other._safe_select_U;
       this->_safe_rank = other._safe_rank;
    }

    ~MyEF(){
       //delete [] _L;
       //delete [] _U;
       //delete [] _M;
    }
};
