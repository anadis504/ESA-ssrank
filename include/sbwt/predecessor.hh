#ifndef _PRED_H_
#define _PRED_H_

#include <stdio.h>
#include <stdlib.h>

#include <chrono>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>
#include <ratio>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;
using namespace std::chrono;

class Predecessor {
   public:
    Predecessor() {
    }
    Predecessor(const vector<uint64_t> &data) {
        _r = data.size();

        // Initialize gaps
        _gaps.reserve(_r);
        _gaps.resize(_r);

        if (_r > 0) {
            // Initialize gaps
            _gaps[0] = data[0];
            uint64_t maxg = 0;
            for (uint64_t i = 1; i < _r; i++) {
                uint64_t g = data[i] - data[i - 1];
                if (g > 65535) {  // bigger than 16 bit integers
                    cerr << "Uh oh: big gap\n";
                }
                if (g > maxg) maxg = g;
                _gaps[i] = (uint16_t)g;
            }
            std::cerr << "max gap: "<<maxg<<'\n';
        }
        _n = 0;
        for (size_t j = 0; j < _r; j++) _n += _gaps[j];
        
        // build predecessor data structure for phrase starting positions
        _numpblocks = (_r / _pbs);
        if (_r % _pbs) {
            _numpblocks++;
        }
        _partialSums.reserve(_numpblocks + 1);
        _partialSums.resize(_numpblocks + 1);
        _partialSums[_numpblocks] = UINT64_MAX;  //_n+1;
        uint64_t sum = 0;
        for (uint64_t i = 0; i < _r; i++) {
            if (i % _pbs == 0) {
                _partialSums[i / _pbs] = sum;
            }
            sum += _gaps[i];
        }
        uint64_t maxdiff = 0;
        uint64_t mindiff = UINT64_MAX;
        for(uint64_t i=1;i<_numpblocks;i++){
           uint64_t d = _partialSums[i] - _partialSums[i-1];
           if(d < mindiff) mindiff = d;
           if(d > maxdiff) maxdiff = d;
        }
        cerr << "mind maxd: "<<mindiff<<' '<<maxdiff<<'\n';

        _msbn = __builtin_clzll(
            _n);  // Returns the number of leading 0-bits in _n,
                  // starting at the most significant bit position.
                  // If _n is 0, the result is undefined.
        _jtShift = (64 - _msbn) - _logjtSize;
        _jumpTable.reserve(_jtSize);
        _jumpTable.resize(_jtSize);

        for (uint32_t j = 0; j < _jtSize; j++) {
            _jumpTable[j] = UINT32_MAX;  // was useful for debugging
        }

        _jumpTable[0] = 0;
        uint32_t jtp = 0;
        uint32_t prevKey = 0;
        uint32_t run = 1, maxRun = 0;
        for (uint32_t i = 1; i < _numpblocks; i++) {
            uint32_t key = getJumpTableKey(_partialSums[i]);
            if (key != prevKey) {
                if (run > maxRun) maxRun = run;
                run = 0;
                while (jtp < key) {
                    _jumpTable[jtp++] = _jumpTable[prevKey];
                }
                _jumpTable[key] = i;
            }
            // jtp++;
            prevKey = key;
            run++;
        }
        while (jtp < _jtSize) {
            _jumpTable[jtp++] = _jumpTable[prevKey];
        }
    }

    Predecessor(Predecessor &other){
       this->_n = other._n;
       this->_r = other._r;

       this->_pbs = other._pbs;
       this->_numpblocks = other._numpblocks;
       this->_partialSums = other._partialSums;
       this->_gaps = other._gaps;

       this->_msbn = other._msbn;
       this->_jtShift = other._jtShift;
       this->_jumpTable = other._jumpTable;
       this->_jtSize = other._jtSize;
       this->_logjtSize = other._logjtSize;
    }

    inline uint32_t getJumpTableKey(uint64_t index) const{
        return (index >> _jtShift);
    }

    // outputs a pair <p,bool>
    //  p is the index of the predecessor in the set
    //  bool is 1 if the value at index p is equal to key, else 0
    pair<int64_t, bool> inline getPred(int64_t key) const{
        if (key < _partialSums[0]) {
            return {-1, 0};
        }
        const uint64_t *x = std::lower_bound(
            _partialSums.data(), _partialSums.data() + _numpblocks,
            key);  // Returns an iterator pointing to the first element in the
                   // range [first, last) that does not satisfy element < value
        if (x - _partialSums.data()) {  // x > _partialSum start
            x--;  // x now pointing at predecessor in _partialSums
        }
        uint64_t sum = *x;
        uint64_t p = _pbs * (x - _partialSums.data());
        while (p < _r) {  // sum < key &&s
            uint64_t len = _gaps[p];
            if (sum + len > key) break;
            p++;
            sum += len;
        }
        return {p - 1, (sum == key)};
    }

    pair<int64_t, bool> inline getPredWithJumpTable(uint64_t pkey) const{
        if (pkey < _partialSums[0]) {
            return {-1, 0};
        }
        if (pkey > _partialSums[_numpblocks]) {
            pkey = _partialSums[_numpblocks];
        }
        uint32_t jkey = getJumpTableKey(pkey);
        const uint64_t *x = _partialSums.data() + _jumpTable[jkey];
        while (*x < pkey) {
            x++;
        }
        if (x - _partialSums.data()) {
            x--;  // x now pointing at predecessor in _partialSums
        }
        uint64_t sum = *x;
        int64_t p =
            _pbs *
            (x -
             _partialSums.data());  // should surely be replaced by a shift (?)
        while (p < _r) {
            uint64_t len = _gaps[p];
            if (sum + len > pkey) {
                break;
            }
            p++;
            sum += len;
        }
        return {p - 1, (sum == pkey)};
    }

    int64_t rank(int64_t pos) const{
       pair<int64_t, bool> r = getPredWithJumpTable(pos);
       //pair<int64_t, bool> r = getPred(pos);
       return (r.first+1)-((uint64_t)(r.second));
    } 

    size_t getn() const { return _n; }

    uint64_t sizeInBytes() const{
        uint64_t sz = 0;
        sz += _r * sizeof(uint16_t);
        sz += _numpblocks * sizeof(uint64_t);
        return sz;
    }

    int64_t serialize(std::ostream& os) const{
        int64_t written = 0;
        os.write((char *)&_n, sizeof(size_t));
        os.write((char *)&_r, sizeof(size_t));
        os.write((char *)&_pbs, sizeof(uint32_t));
        os.write((char *)&_numpblocks, sizeof(uint32_t));
        os.write((char *)&_msbn, sizeof(uint32_t));
        os.write((char *)&_jtShift, sizeof(uint32_t));
        os.write((char *)&_jtSize, sizeof(uint32_t));
        os.write((char *)&_logjtSize, sizeof(uint32_t));
        written += sizeof(uint32_t) * 6 + sizeof(size_t) * 2;

        //now serialize vectors _partialSums, _gaps, _jumpTable
        
        size_t size = _partialSums.size();
       std::cerr << "size: "<<size<<'\n';
        os.write(reinterpret_cast<char const*>(&size), sizeof(size));
        os.write(reinterpret_cast<char const*>(_partialSums.data()), _partialSums.size() * sizeof(uint64_t));
        written += sizeof(size) + _partialSums.size() * sizeof(uint64_t);

        size = _gaps.size();
       std::cerr << "size: "<<size<<'\n';
        os.write(reinterpret_cast<char const*>(&size), sizeof(size));
        os.write(reinterpret_cast<char const*>(_gaps.data()), _gaps.size() * sizeof(uint16_t));
        written += sizeof(size) + _gaps.size() * sizeof(uint16_t);

        size = _jumpTable.size();
       std::cerr << "size: "<<size<<'\n';
        os.write(reinterpret_cast<char const*>(&size), sizeof(size));
        os.write(reinterpret_cast<char const*>(_jumpTable.data()), _jumpTable.size() * sizeof(uint32_t));
        written += sizeof(size) + _jumpTable.size() * sizeof(uint32_t);
        
        return written;
    }

    void load(std::istream& is){
       is.read((char *)&_n, sizeof(size_t));
       is.read((char *)&_r, sizeof(size_t));
       is.read((char *)&_pbs, sizeof(uint32_t));
       is.read((char *)&_numpblocks, sizeof(uint32_t));
       is.read((char *)&_msbn, sizeof(uint32_t));
       is.read((char *)&_jtShift, sizeof(uint32_t));
       is.read((char *)&_jtSize, sizeof(uint32_t));
       is.read((char *)&_logjtSize, sizeof(uint32_t));

       std::cerr << "_n _r _pbs _numpblocks _msbn _jtShift _jtSize _logjtSize: "<<_n<<' '<<_r<<' '<<_pbs<<' '<<_numpblocks<<' '<<_msbn<<' '<<_jtShift<<' '<<_jtSize<<' '<<_logjtSize<<'\n';

       //now load data for vectors _partialSums, _gaps, _jumpTable
       size_t size;
       is.read(reinterpret_cast<char*>(&size), sizeof(size));
       std::cerr << "size: "<<size<<'\n';
       _partialSums.reserve(size);
       _partialSums.resize(size);
       is.read(reinterpret_cast<char*>(_partialSums.data()), _partialSums.size() * sizeof(uint64_t));

       is.read(reinterpret_cast<char*>(&size), sizeof(size));
       std::cerr << "size: "<<size<<'\n';
       _gaps.resize(size);
       is.read(reinterpret_cast<char*>(_gaps.data()), _gaps.size() * sizeof(uint16_t));

       is.read(reinterpret_cast<char*>(&size), sizeof(size));
       std::cerr << "size: "<<size<<'\n';
       _jumpTable.resize(size);
       is.read(reinterpret_cast<char*>(_jumpTable.data()), _jumpTable.size() * sizeof(uint32_t));
    }    

   private:
    inline uint64_t getGap(uint64_t index) { return (_gaps[index]); }

    size_t _n = 0;  // length of BWT and text
    size_t _r = 0;  // number of runs

    uint32_t _pbs = 64;             // predecessor block size
    uint32_t _numpblocks = 0;       // number of blocks in predecessor structure
    vector<uint64_t> _partialSums;  // every _pbs^th
    vector<uint16_t> _gaps;

    uint32_t _msbn;
    uint32_t _jtShift;
    vector<uint32_t> _jumpTable;
    uint32_t _jtSize = 131072;
    uint32_t _logjtSize = 17;
    //uint32_t _jtSize = 65536;
    //uint32_t _logjtSize = 16;
};

#endif
