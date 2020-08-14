// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#pragma once
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>
#include "blockRadixSort.h"
#include "graph.h"
#include "parallel.h"
#include "quickSort.h"
#include "utils.h"
//#include <pmem_allocator.h>
#include <libpmem.h>
#include "nodePartitioner.h"
using namespace std;

typedef pair<uintE, uintE> intPair;
typedef pair<uintE, pair<uintE, intE> > intTriple;

template <class E>
struct pairFirstCmp {
  bool operator()(pair<uintE, E> a, pair<uintE, E> b) {
    return a.first < b.first;
  }
};

template <class E>
struct getFirst {
  uintE operator()(pair<uintE, E> a) { return a.first; }
};

template <class IntType>
struct pairBothCmp {
  bool operator()(pair<uintE, IntType> a, pair<uintE, IntType> b) {
    if (a.first != b.first) return a.first < b.first;
    return a.second < b.second;
  }
};

/*
// A structure that keeps a sequence of strings all allocated from
// the same block of memory
struct words {
  long n; // total number of characters
  char* Chars;  // array storing all strings
  long m; // number of substrings
  char** Strings; // pointers to strings (all should be null terminated)
  words() {}
words(char* C, long nn, char** S, long mm)
: Chars(C), n(nn), Strings(S), m(mm) {}
  void del() {free(Chars); free(Strings);}
};
*/

inline bool isSpace(char c) {
  switch (c) {
    case '\r':
    case '\t':
    case '\n':
    case 0:
    case ' ':
      return true;
    default:
      return false;
  }
}

_seq<char> mmapStringFromFile(const char* filename) {
  struct stat sb;
  int fd = open(filename, O_RDONLY);
  if (fd == -1) {
    perror("open");
    exit(-1);
  }
  if (fstat(fd, &sb) == -1) {
    perror("fstat");
    exit(-1);
  }
  if (!S_ISREG(sb.st_mode)) {
    perror("not a file\n");
    exit(-1);
  }
  char* p =
      static_cast<char*>(mmap(0, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0));
  if (p == MAP_FAILED) {
    perror("mmap");
    exit(-1);
  }
  if (close(fd) == -1) {
    perror("close");
    exit(-1);
  }
  size_t n = sb.st_size;
  //  char *bytes = newA(char, n);
  //  parallel_for(size_t i=0; i<n; i++) {
  //    bytes[i] = p[i];
  //  }
  //  if (munmap(p, sb.st_size) == -1) {
  //    perror("munmap");
  //    exit(-1);
  //  }
  //  cout << "mmapped" << endl;
  //  free(bytes);
  //  exit(0);
  return _seq<char>(p, n);
}

_seq<char> readStringFromFile(char* fileName) {
  ifstream file(fileName, ios::in | ios::binary | ios::ate);
  if (!file.is_open()) {
    std::cout << "Unable to open file: " << fileName << std::endl;
    abort();
  }
  long end = file.tellg();
  file.seekg(0, ios::beg);
  long n = end - file.tellg();
  char* bytes = newA(char, n + 1);
  file.read(bytes, n);
  file.close();
  return _seq<char>(bytes, n);
}

// parallel code for converting a string to words
words stringToWords(char* Str, long n) {
  { parallel_for(long i = 0; i < n; i++) if (isSpace(Str[i])) Str[i] = 0; }

  // mark start of words
  bool* FL = newA(bool, n);
  FL[0] = Str[0];
  { parallel_for(long i = 1; i < n; i++) FL[i] = Str[i] && !Str[i - 1]; }

  // offset for each start of word
  _seq<long> Off = sequence::packIndex<long>(FL, n);
  long m = Off.n;
  long* offsets = Off.A;

  // pointer to each start of word
  char** SA = newA(char*, m);
  { parallel_for(long j = 0; j < m; j++) SA[j] = Str + offsets[j]; }

  free(offsets);
  free(FL);
  return words(Str, n, SA, m);
}

// re write this to get the largest idx possible
//------------
size_t binarySearch(uintT* a, size_t low, size_t high, size_t m) {
  if ((high - low) <= 1) {
    if ((a[low + 1] > m) && (a[low] <= m)) {
      return low;
    } else {
      std::cerr << "binsearch error: "
                << " l " << low << " h " << high << " m " << m << "\n";
      exit(1);
    }
  } else if (high > low) {
    size_t mid = low + (high - low) / 2;
    if (a[mid] > m) return binarySearch(a, low, mid, m);
    return binarySearch(a, mid, high, m);
  } else {
    std::cerr << "binsearch error: "
              << " l " << low << " h " << high << " m " << m << "\n";
    exit(1);
  }
}
//------------

template <class vertex>
graph<vertex> readGraphFromFile(char* fname, bool isSymmetric, bool mmap,
                                long capacity = 0, char* npf = nullptr,
                                bool rmetis = false) {
  words W;
  if (mmap) {
    _seq<char> S = mmapStringFromFile(fname);
    char* bytes = newA(char, S.n);
    // Cannot mutate the graph unless we copy.
    parallel_for(size_t i = 0; i < S.n; i++) { bytes[i] = S.A[i]; }
    if (munmap(S.A, S.n) == -1) {
      perror("munmap");
      exit(-1);
    }
    S.A = bytes;
    W = stringToWords(S.A, S.n);
  } else {
    _seq<char> S = readStringFromFile(fname);
    W = stringToWords(S.A, S.n);
  }
#ifndef WEIGHTED
  if (W.Strings[0] != (string) "AdjacencyGraph") {
#else
  if (W.Strings[0] != (string) "WeightedAdjacencyGraph") {
#endif
    cout << "Bad input file l0" << endl;
    abort();
  }

  long len = W.m - 1;
  long n = atol(W.Strings[1]);
  long m = atol(W.Strings[2]);
#ifndef WEIGHTED
  if (len != n + m + 2) {
#else
  if (len != n + 2 * m + 2) {
#endif
    cout << "Bad input file count" << endl;
    abort();
  }

  string path = "/mnt/pmem0/mina";

  char* pmemaddr;
  int is_pmem;
  size_t mapped_len;
  ifstream f(path.c_str());

  const static size_t PMEMGRAPHSIZE = 200000000000;

  if (f.good()) {
    pmemaddr =
        (char*)pmem_map_file(path.c_str(), PMEMGRAPHSIZE, PMEM_FILE_CREATE,
                             0666, &mapped_len, &is_pmem);
  } else {
    pmemaddr = (char*)pmem_map_file(path.c_str(), PMEMGRAPHSIZE,
                                    PMEM_FILE_CREATE | PMEM_FILE_EXCL, 0666,
                                    &mapped_len, &is_pmem);
  }

  if (pmemaddr == NULL) {
    perror("pmem_map_file");
    exit(1);
  }

  if (is_pmem) {
    cout << "Pmem successfully allocated" << endl;
  } else {
    cout << "Error not pmem" << endl;
    exit(1);
  }

  uintE* nv = reinterpret_cast<uintE*>(pmemaddr);

  // libmemkind::pmem::allocator<uintE> alc{ "/mnt/pmem12", 102400 };
  // std::allocator<uintE> alc;

  // increase size of offsets to n+1
  uintT* offsets = newA(uintT, n + 1);

#ifndef WEIGHTED
  uintE* edges = newA(uintE, m);

  // memmove(&nv[0], pmemaddr, sizeof(nv))
  // uintE* nv = alc.allocate(m);
#else
  intE* edges = newA(intE, 2 * m);
#endif

  {
    parallel_for(long i = 0; i < n; i++) {
      offsets[i] = atol(W.Strings[i + 3]);
      // if ((i == 37296) || (i == 37302) || (i== 37303))
      //  std::cout << "offsets og " << offsets[i] << " and atol from " << i+3
      //  << std::endl;
    }
  }
  // can remove the if conditions
  offsets[n] = m;

  NodePartitioner* np;
  uintT* tOffsets = newA(uintE, n);
  memset(tOffsets, 0, n * sizeof(uintT));
  if (npf != nullptr) {
    _seq<char> dn = readStringFromFile(npf);
    words dnw = stringToWords(dn.A, dn.n);
    // std::cerr << "got to this point: dn | " << dn.n << " " << (bool) dn.A[0]
    // << std::endl;

    // for (size_t i = 0; i<n; i++)
    //  std::cout << "NP| " << i << " " << np->dram_nodes[i] << " " <<
    //  np->if_dram(i) << " : " << np->partition_offset(i) << std::endl;
    if (!isSymmetric) {
      parallel_for(int i = 0; i < m; i++) {
        writeAdd(&tOffsets[atol(W.Strings[i + n + 3])], (uintT)1);
      }

      np = new NodePartitioner(n, m, offsets, tOffsets, dnw, rmetis);
      long collect_sum = 0;
      for (int i = 0; i < n; i++) {
        long tmp = collect_sum;
        collect_sum += tOffsets[i];
        tOffsets[i] = tmp;
      }
    } else {
      np = new NodePartitioner(n, m, offsets, dnw, rmetis);
    }
  } else {
    if (!isSymmetric) {
      parallel_for(int i = 0; i < m; i++) {
        writeAdd(&tOffsets[atol(W.Strings[i + n + 3])], (uintT)1);
      }

      np = new NodePartitioner(n, m, offsets, tOffsets, capacity);
      long collect_sum = 0;
      for (int i = 0; i < n; i++) {
        long tmp = collect_sum;
        collect_sum += tOffsets[i];
        tOffsets[i] = tmp;
      }
    } else {
      np = new NodePartitioner(n, m, offsets, capacity);
    }
  }

  // for debugging
  size_t edge_count = 0;
  size_t node_count = 0;
  for (size_t i = 0; i < n; i++) {
    if (np->if_dram(i)) {
      edge_count += ((i == n - 1) ? m : offsets[i + 1]) - offsets[i];
      node_count++;
    }
  }
  if (isSymmetric) {
    std::cout << "Percent of all nodes on dram: " << node_count << "/" << n
              << std::endl;
    std::cout << "Percent of all edges on dram: " << edge_count << "/" << m
              << std::endl;
  } else {
    std::cout << "Percent of all (out) nodes on dram: " << node_count << "/"
              << n << std::endl;
    std::cout << "Percent of all (out) edges on dram: " << edge_count << "/"
              << m << std::endl;
    size_t in_edge_count = 0;
    size_t in_node_count = 0;
    for (size_t i = 0; i < n; i++) {
      if (np->if_dram(i + n)) {
        in_edge_count += ((i == n - 1) ? m : tOffsets[i + 1]) - tOffsets[i];
        in_node_count++;
      }
    }
    std::cout << "Percent of all (in) nodes on dram: " << in_node_count << "/"
              << n << std::endl;
    std::cout << "Percent of all (in) edges on dram: " << in_edge_count << "/"
              << m << std::endl;
  }

  // for (size_t i = 0; i<n; i++)
  //  std::cout << "NP| " << i << " " << np->if_dram(i) << " : " <<
  //  np->partition_offset(i) << std::endl; why not marked as volatile???
  uintT dram_partition_split = 0;

  {
    parallel_for(long i = 0; i < m; i++) {
#ifndef WEIGHTED
      // size_t binary_search(Sequence I, typename Sequence::T v, const F& less)
      // {
      size_t j = binarySearch(offsets, 0, n, i);  // rename this better
      // while(i > offsets[j]) j++;
      uintT po = np->partition_offset(j);
      uintT l = ((j == n - 1) ? m : offsets[j + 1]) - offsets[j];
      // if ((i == 112442) || (j == 37302) || (i== 112441) || (i== 112440))
      //  std::cout << " n " << j << " m " << i << " l " << l << " offsets " <<
      //  offsets[j] << " " << offsets[j+1] << std::endl;
      if (np->if_dram(j)) {
        edges[po + (i - offsets[j])] = atol(W.Strings[i + n + 3]);
        // if ((i == 112442) || (j == 37302) || (i== 112441))
        //  std::cerr << "writing edges | n " << j << " m:" << i << " d " <<
        //  W.Strings[i+n+3] \
      //    << " from idx " << i+n+3 << " poff: " << np->partition_offset(j) +
        //    (i - offsets[j])<< std::endl;
        pbbs::write_min(&dram_partition_split, po + l, std::greater<uintT>());
        // std::cout << "dram partition split new | " << dram_partition_split <<
        // " = " << po << "+" << l << std::endl;
      } else {
        // if ((i == 112442) || (j == 37302) || (i== 112441))
        //  std::cerr << "writing nv edges |n " << j << " m:" << i << " d " <<
        //  W.Strings[i+n+3] \
      //  << " poff: " << np->partition_offset(j) + (i - offsets[j]) <<
        //  std::endl;
        nv[po + (i - offsets[j])] = atol(W.Strings[i + n + 3]);
      }
#else
      edges[2 * i] = atol(W.Strings[i + n + 3]);
      edges[2 * i + 1] = atol(W.Strings[i + n + m + 3]);
#endif
    }
  }
  // W.del(); // to deal with performance bug in malloc

  // for(size_t i=0;i<2*m;i++) std::cout << "edges check " << edges[i] << " a:"
  // << &edges[i] << std::endl; for(size_t i=0;i<2*m;i++) std::cout << "nv edges
  // check " << nv[i] << " a:" << &nv[i] << std::endl; std::cout << "Made
  // edges\n";

  vertex* v = newA(vertex, n);

  {
    parallel_for(uintT i = 0; i < n; i++) {
      uintT o = offsets[i];
      uintT l = ((i == n - 1) ? m : offsets[i + 1]) - offsets[i];
      v[i].setOutDegree(l);
#ifndef WEIGHTED
      if (np->if_dram(i)) {
        v[i].setOutNeighbors(edges + np->partition_offset(i));
        // TODO fix the offsets or array of pointer problem
        // if (pmem_is_pmem(edges+o, l*sizeof(uintE))) { std::cerr << "pmem when
        // its not" << i << " " << (void*)(edges+o) << std::endl; } if (i ==
        // 37302)
        //  std::cout << "added out neighbor dram| " << i << " " <<
        //  np->partition_offset(i) << " with degr " << l << \
      //    " from " << offsets[i+1] << " - " << offsets[i] << std::endl;
        // if (i<5) std::cout << "first outNeighbor dram| " << i << " " <<
        // v[i].getOutNeighbor(0) << std::endl;
      } else {
        v[i].setOutNeighbors(nv + np->partition_offset(i));
        // if (!pmem_is_pmem(nv+o, l*sizeof(uintE))) { std::cerr << "not pmem "
        // << i << " " << (void*)(nv+o) << std::endl; } std::cout << "added out
        // neighbor nvram| " << i << " " << np->partition_offset(i) <<
        // std::endl; if (i<5) std::cout << "first outNeighbor nvram| " << i <<
        // " " << v[i].getOutNeighbor(0) << std::endl;
      }
#else
      v[i].setOutNeighbors(edges + 2 * o);
#endif
    }
  }
  // std::cout << "Made outneighbors\n";

  // for(size_t i=0;i<2*m;i++) std::cout << "edges check " << edges[i] << " a:"
  // << &edges[i] << std::endl; for(size_t i=0;i<2*m;i++) std::cout << "nv edges
  // check " << nv[i] << " a:" << &nv[i] << std::endl; for (size_t i=0; i<n;
  // i++)
  // {
  //  for (size_t j=0; j<v[i].getOutDegree(); j++)
  //    std::cout << "outneighbors check " << i << "-" << j << " : " <<
  //    v[i].getOutNeighbor(j) << " a:" << v[i].getOutNeighbor(j) << std::endl;
  //}

  if (!isSymmetric) {
    // this should all be moved to nvram
    // rn it doesnt matter that its on dram but it should be nvram
    // uintT* tOffsets = newA(uintT,n);
    //{parallel_for(long i=0;i<n;i++) tOffsets[i] = INT_T_MAX;}
#ifndef WEIGHTED
    intPair* temp = newA(intPair, m);
#else
    intTriple* temp = newA(intTriple, m);
#endif
    {
      parallel_for(long i = 0; i < n; i++) {
        uintT o = offsets[i];
        for (uintT j = 0; j < v[i].getOutDegree(); j++) {
#ifndef WEIGHTED
          temp[o + j] = make_pair(v[i].getOutNeighbor(j), i);
          // std::cout << "temp init " << j << " : " << v[i].getOutNeighbor(j)
          // << ", " << i << std::endl;
          assert(v[i].getOutNeighbor(j) < 16777216);
#else
          temp[o + j] = make_pair(v[i].getOutNeighbor(j),
                                  make_pair(i, v[i].getOutWeight(j)));
#endif
        }
      }
    }
    free(offsets);

#ifndef WEIGHTED
#ifndef LOWMEM
    intSort::iSort(temp, m, n + 1, getFirst<uintE>());
#else
    quickSort(temp, m, pairFirstCmp<uintE>());
#endif
#else
#ifndef LOWMEM
    intSort::iSort(temp, m, n + 1, getFirst<intPair>());
#else
    quickSort(temp, m, pairFirstCmp<intPair>());
#endif
#endif

    std::cout << "Made sort\n";

    // tOffsets[temp[0].first] = 0;
#ifndef WEIGHTED
    uintE* inEdges = newA(uintE, m);
    // inEdges[0] = temp[0].second;

    // uintE* niv = alc.allocate(m);

    // uintE* niv = bufaddr + m + 2;
    // niv[0] = temp[0].second;
    // for(size_t i=0;i<m;i++) std::cout << "edges check " << inEdges[i] << "
    // a:" << &(inEdges[i]) << std::endl;
#else
    intE* inEdges = newA(intE, 2 * m);
    inEdges[0] = temp[0].second.first;
    inEdges[1] = temp[0].second.second;
#endif
    // for (size_t i=0; i<m; i++) std::cout << "temp " << temp[i].first << " "
    // << temp[i].second << std::endl; std::cout << "dram_partition_split " <<
    // dram_partition_split << std::endl;
    {
      parallel_for(long i = 0; i < m; i++) {  //
#ifndef WEIGHTED
        size_t j = temp[i].first;
        if (np->if_dram(j + n)) {
          // std::cerr << "writing edges | n " << j << " m:" << i << " ";
          // std::cerr << "d " << temp[i].second << " poff: " <<
          // np->partition_offset(j+n) + \
        //  (i - tOffsets[j]) - dram_partition_split << std::endl;
          // std::cerr << " po: " << np->partition_offset(j+n) << " i: " << i <<
          // "tOffsets[j]: " << \
        //  tOffsets[j] << " dps: " << dram_partition_split << std::endl;
          inEdges[np->partition_offset(j + n) + (i - tOffsets[j]) -
                  dram_partition_split] = temp[i].second;
        } else {
          // std::cerr << "writing nv edges |n " << j << " m:" << i <<  " ";
          // std::cerr << " d " << temp[i].second << " poff: " <<
          // np->partition_offset(j+n) + \
        //  (i - tOffsets[j]) << " po: " << np->partition_offset(j+n) << " i:
          //  " << i << "tOffsets[j]: " << tOffsets[j] << std::endl;
          nv[np->partition_offset(j + n) + (i - tOffsets[j])] = temp[i].second;
        }
        // inEdges[i] = temp[i].second;
        // niv[i] = temp[i].second;
#else
        inEdges[2 * i] = temp[i].second.first;
        inEdges[2 * i + 1] = temp[i].second.second;
#endif
        // if(temp[i].first != temp[i-1].first) {
        //  tOffsets[temp[i].first] = i;
        //}
      }
    }
    // std::cout << "Made inedges\n";
    free(temp);

    // fill in offsets of degree 0 vertices by taking closest non-zero
    // offset to the right
    // sequence::scanIBack(tOffsets,tOffsets,n,minF<uintT>(),(uintT)m);

    // why cant I run this with one thread?
    {
      parallel_for(long i = 0; i < n; i++) {
        uintT o = tOffsets[i];
        uintT l = ((i == n - 1) ? m : tOffsets[i + 1]) - tOffsets[i];
        // std::cout << "tOffsets " << i << ":" << l << std::endl;
        v[i].setInDegree(l);
#ifndef WEIGHTED
        // v[i].setInNeighbors(inEdges+o);
        // v[i].setInNeighbors(niv+o);
        if (np->if_dram(i + n)) {
          v[i].setInNeighbors(inEdges + np->partition_offset(i + n) -
                              dram_partition_split);
          // TODO fix the offsets or array of pointer problem
          // if (pmem_is_pmem(edges+o, l*sizeof(uintE))) { std::cerr << "pmem
          // when its not" << i << " " << (void*)(edges+o) << std::endl; }
          // std::cout << "added in neighbor dram| " << i << " " <<
          // np->partition_offset(i+n)-dram_partition_split << std::endl; if
          // (i<5) std::cout << "first outNeighbor dram| " << i << " " <<
          // v[i].getOutNeighbor(0) << std::endl;
        } else {
          v[i].setInNeighbors(nv + np->partition_offset(i + n));
          // if (!pmem_is_pmem(nv+o, l*sizeof(uintE))) { std::cerr << "not pmem
          // " << i << " " << (void*)(nv+o) << std::endl; } std::cout << "added
          // in neighbor nvram| " << i << " " << np->partition_offset(i+n) <<
          // std::endl; if (i<5) std::cout << "first outNeighbor nvram| " << i
          // << " " << v[i].getOutNeighbor(0) << std::endl;
        }
#else
        v[i].setInNeighbors(inEdges + 2 * o);
#endif
      }
    }
    // for(size_t i=0;i<m;i++) std::cout << "edges check " << inEdges[i] << "
    // a:" << &(inEdges[i]) << std::endl; for(size_t i=0;i<2*m;i++) std::cout <<
    // "nv edges check " << nv[i] << " a:" << &nv[i] << std::endl; for (size_t
    // i=0; i<n; i++) {
    //  for (size_t j=0; j<v[i].getInDegree(); j++)
    //    std::cout << "inneighbors check " << i << "-" << j << " : " <<
    //    v[i].getInNeighbor(j) << std::endl;
    //}
    // std::cout << "Made inneighbors\n";
    free(tOffsets);
    Uncompressed_Mem<vertex>* mem =
        new Uncompressed_Mem<vertex>(v, n, m, edges, inEdges);
    return graph<vertex>(v, n, m, mem);
  } else {
    free(offsets);
    Uncompressed_Mem<vertex>* mem =
        new Uncompressed_Mem<vertex>(v, n, m, edges);
    return graph<vertex>(v, n, m, mem);
  }
}

template <class vertex>
graph<vertex> readGraphFromBinary(char* iFile, bool isSymmetric) {
  char* config = (char*)".config";
  char* adj = (char*)".adj";
  char* idx = (char*)".idx";
  char configFile[strlen(iFile) + strlen(config) + 1];
  char adjFile[strlen(iFile) + strlen(adj) + 1];
  char idxFile[strlen(iFile) + strlen(idx) + 1];
  *configFile = *adjFile = *idxFile = '\0';
  strcat(configFile, iFile);
  strcat(adjFile, iFile);
  strcat(idxFile, iFile);
  strcat(configFile, config);
  strcat(adjFile, adj);
  strcat(idxFile, idx);

  ifstream in(configFile, ifstream::in);
  long n;
  in >> n;
  in.close();

  ifstream in2(adjFile, ifstream::in | ios::binary);  // stored as uints
  in2.seekg(0, ios::end);
  long size = in2.tellg();
  in2.seekg(0);
#ifdef WEIGHTED
  long m = size / (2 * sizeof(uint));
#else
  long m = size / sizeof(uint);
#endif
  char* s = (char*)malloc(size);
  in2.read(s, size);
  in2.close();
  uintE* edges = (uintE*)s;

  ifstream in3(idxFile, ifstream::in | ios::binary);  // stored as longs
  in3.seekg(0, ios::end);
  size = in3.tellg();
  in3.seekg(0);
  if (n != size / sizeof(intT)) {
    cout << "File size wrong\n";
    abort();
  }

  char* t = (char*)malloc(size);
  in3.read(t, size);
  in3.close();
  uintT* offsets = (uintT*)t;

  vertex* v = newA(vertex, n);
#ifdef WEIGHTED
  intE* edgesAndWeights = newA(intE, 2 * m);
  {
    parallel_for(long i = 0; i < m; i++) {
      edgesAndWeights[2 * i] = edges[i];
      edgesAndWeights[2 * i + 1] = edges[i + m];
    }
  }
  // free(edges);
#endif
  {
    parallel_for(long i = 0; i < n; i++) {
      uintT o = offsets[i];
      uintT l = ((i == n - 1) ? m : offsets[i + 1]) - offsets[i];
      v[i].setOutDegree(l);
#ifndef WEIGHTED
      v[i].setOutNeighbors((uintE*)edges + o);
#else
      v[i].setOutNeighbors(edgesAndWeights + 2 * o);
#endif
    }
  }

  if (!isSymmetric) {
    uintT* tOffsets = newA(uintT, n);
    { parallel_for(long i = 0; i < n; i++) tOffsets[i] = INT_T_MAX; }
#ifndef WEIGHTED
    intPair* temp = newA(intPair, m);
#else
    intTriple* temp = newA(intTriple, m);
#endif
    {
      parallel_for(intT i = 0; i < n; i++) {
        uintT o = offsets[i];
        for (uintT j = 0; j < v[i].getOutDegree(); j++) {
#ifndef WEIGHTED
          temp[o + j] = make_pair(v[i].getOutNeighbor(j), i);
#else
          temp[o + j] = make_pair(v[i].getOutNeighbor(j),
                                  make_pair(i, v[i].getOutWeight(j)));
#endif
        }
      }
    }
    free(offsets);
#ifndef WEIGHTED
#ifndef LOWMEM
    intSort::iSort(temp, m, n + 1, getFirst<uintE>());
#else
    quickSort(temp, m, pairFirstCmp<uintE>());
#endif
#else
#ifndef LOWMEM
    intSort::iSort(temp, m, n + 1, getFirst<intPair>());
#else
    quickSort(temp, m, pairFirstCmp<intPair>());
#endif
#endif
    tOffsets[temp[0].first] = 0;
#ifndef WEIGHTED
    uintE* inEdges = newA(uintE, m);
    inEdges[0] = temp[0].second;
#else
    intE* inEdges = newA(intE, 2 * m);
    inEdges[0] = temp[0].second.first;
    inEdges[1] = temp[0].second.second;
#endif
    {
      parallel_for(long i = 1; i < m; i++) {
#ifndef WEIGHTED
        inEdges[i] = temp[i].second;
#else
        inEdges[2 * i] = temp[i].second.first;
        inEdges[2 * i + 1] = temp[i].second.second;
#endif
        if (temp[i].first != temp[i - 1].first) {
          tOffsets[temp[i].first] = i;
        }
      }
    }
    free(temp);
    // fill in offsets of degree 0 vertices by taking closest non-zero
    // offset to the right
    sequence::scanIBack(tOffsets, tOffsets, n, minF<uintT>(), (uintT)m);
    {
      parallel_for(long i = 0; i < n; i++) {
        uintT o = tOffsets[i];
        uintT l = ((i == n - 1) ? m : tOffsets[i + 1]) - tOffsets[i];
        v[i].setInDegree(l);
#ifndef WEIGHTED
        v[i].setInNeighbors((uintE*)inEdges + o);
#else
        v[i].setInNeighbors((intE*)(inEdges + 2 * o));
#endif
      }
    }
    free(tOffsets);
#ifndef WEIGHTED
    Uncompressed_Mem<vertex>* mem =
        new Uncompressed_Mem<vertex>(v, n, m, edges, inEdges);
    return graph<vertex>(v, n, m, mem);
#else
    Uncompressed_Mem<vertex>* mem =
        new Uncompressed_Mem<vertex>(v, n, m, edgesAndWeights, inEdges);
    return graph<vertex>(v, n, m, mem);
#endif
  }
  free(offsets);
#ifndef WEIGHTED
  Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(v, n, m, edges);
  return graph<vertex>(v, n, m, mem);
#else
  Uncompressed_Mem<vertex>* mem =
      new Uncompressed_Mem<vertex>(v, n, m, edgesAndWeights);
  return graph<vertex>(v, n, m, mem);
#endif
}

template <class vertex>
graph<vertex> readGraph(char* iFile, bool compressed, bool symmetric,
                        bool binary, bool mmap, long capacity = 0,
                        char* npf = nullptr, bool rmetis = false) {
  if (binary)
    return readGraphFromBinary<vertex>(iFile, symmetric);
  else
    return readGraphFromFile<vertex>(iFile, symmetric, mmap, capacity, npf,
                                     rmetis);
}

template <class vertex>
graph<vertex> readCompressedGraph(char* fname, bool isSymmetric, bool mmap) {
  char* s;
  if (mmap) {
    _seq<char> S = mmapStringFromFile(fname);
    // Cannot mutate graph unless we copy.
    char* bytes = newA(char, S.n);
    parallel_for(size_t i = 0; i < S.n; i++) { bytes[i] = S.A[i]; }
    if (munmap(S.A, S.n) == -1) {
      perror("munmap");
      exit(-1);
    }
    s = bytes;
  } else {
    ifstream in(fname, ifstream::in | ios::binary);
    in.seekg(0, ios::end);
    long size = in.tellg();
    in.seekg(0);
    cout << "size = " << size << endl;
    s = (char*)malloc(size);
    in.read(s, size);
    in.close();
  }

  long* sizes = (long*)s;
  long n = sizes[0], m = sizes[1], totalSpace = sizes[2];

  cout << "n = " << n << " m = " << m << " totalSpace = " << totalSpace << endl;
  cout << "reading file..." << endl;

  uintT* offsets = (uintT*)(s + 3 * sizeof(long));
  long skip = 3 * sizeof(long) + (n + 1) * sizeof(intT);
  uintE* Degrees = (uintE*)(s + skip);
  skip += n * sizeof(intE);
  uchar* edges = (uchar*)(s + skip);

  uintT* inOffsets;
  uchar* inEdges;
  uintE* inDegrees;
  if (!isSymmetric) {
    skip += totalSpace;
    uchar* inData = (uchar*)(s + skip);
    sizes = (long*)inData;
    long inTotalSpace = sizes[0];
    cout << "inTotalSpace = " << inTotalSpace << endl;
    skip += sizeof(long);
    inOffsets = (uintT*)(s + skip);
    skip += (n + 1) * sizeof(uintT);
    inDegrees = (uintE*)(s + skip);
    skip += n * sizeof(uintE);
    inEdges = (uchar*)(s + skip);
  } else {
    inOffsets = offsets;
    inEdges = edges;
    inDegrees = Degrees;
  }

  vertex* V = newA(vertex, n);
  parallel_for(long i = 0; i < n; i++) {
    long o = offsets[i];
    uintT d = Degrees[i];
    V[i].setOutDegree(d);
    V[i].setOutNeighbors(edges + o);
  }

  if (sizeof(vertex) == sizeof(compressedAsymmetricVertex)) {
    parallel_for(long i = 0; i < n; i++) {
      long o = inOffsets[i];
      uintT d = inDegrees[i];
      V[i].setInDegree(d);
      V[i].setInNeighbors(inEdges + o);
    }
  }

  cout << "creating graph..." << endl;
  Compressed_Mem<vertex>* mem = new Compressed_Mem<vertex>(V, s);

  graph<vertex> G(V, n, m, mem);
  return G;
}
