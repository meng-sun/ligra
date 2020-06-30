#ifndef NODE_PARTITIONER
#define NODE_PARTITIONER
#include "words.h"
#include "parallel.h"
#include <cassert>
class NodePartitioner{

  uintE* partition_offsets;
  long _n;

public:
  bool* dram_nodes;
  NodePartitioner(): dram_nodes(nullptr), partition_offsets(nullptr), _n(0) {};
  NodePartitioner(const NodePartitioner& np);
  NodePartitioner(size_t n, size_t m, uintE* o, size_t capacity);
  NodePartitioner(size_t n, size_t m, uintE* o, uintT* indo, size_t capacity);
  NodePartitioner(size_t n, size_t m, uintE* o, words& dn, bool rmetis);
  NodePartitioner(size_t n, size_t m, uintE* o,\
     uintT* indo, words& dn, bool rmetis);
  ~NodePartitioner();
  const bool if_dram(size_t i) { return dram_nodes[i]; }
  const uintE partition_offset(size_t i) { return partition_offsets[i]; }

};

NodePartitioner::NodePartitioner(const NodePartitioner& np) {
  _n = np._n;
  dram_nodes = newA(bool,_n);
  partition_offsets = newA(uintE, _n);
  memcpy(dram_nodes, np.dram_nodes, _n);
  //std::cout << "copying NP"<< std::endl;
  //for (size_t i=0; i<_n; i++)
  //  std::cout <<"dram node pre |" << dram_nodes[i] << std::endl;
  memcpy(partition_offsets, np.partition_offsets, _n*sizeof(uintE));
}

NodePartitioner::NodePartitioner(size_t n, size_t m, uintE* o,\
    words& dn, bool rmetis=false) {
  _n = n;

  dram_nodes = newA(bool, n);
  memset(dram_nodes, false, n);

  assert(dn.m == n);
  std::cout << "Using metis partitioning...\n";
  if (rmetis) {
    parallel_for(long i=0; i<n; i++) {
      if (dn.Strings[i] != (string) "1") dram_nodes[i] = true;
    }
  } else {
    parallel_for(long i=0; i<n; i++) {
      //std::cout << "i " << i << ": " << dn.Strings[i] << "\n";
      if (dn.Strings[i] == (string) "1") dram_nodes[i] = true;
      //std::cout << "dram i " << i << ": " << dram_nodes[i] << "\n";
    }
  }

  std::cout << "Counting partition offsets...\n";
  partition_offsets = newA(uintE, n);
  uintE dram_count = 0;
  uintE other_count = 0;
  for (size_t i=0; i<n; i++) {
    if (dram_nodes[i]) {
      partition_offsets[i] = dram_count;
      dram_count += ((i == n-1) ? m : o[i+1])-o[i];
    } else {
      partition_offsets[i] = other_count;
      other_count += ((i == n-1) ? m : o[i+1])-o[i];
    }
  }

  //for (size_t i=0; i<_n; i++)
  //  std::cout <<"dram node pre |" << dram_nodes[i] << std::endl;
  std::cout << "created node partitioner from file" << std::endl;
}

NodePartitioner::NodePartitioner(size_t n, size_t m, uintE* o,\
     uintT* indo, words& dn, bool rmetis=false) {
  _n = n;
  dram_nodes = newA(bool, 2*n);
  memset(dram_nodes, false, 2*n);
  partition_offsets = newA(uintE, 2*n);

  std::cout << "Using metis partitioning...\n";
  std::cout << "Using asymmetric graph..." << std::endl;

  assert(dn.m == n);

  if (rmetis) {
    parallel_for(long i=0; i<n; i++) {
      if (dn.Strings[i] != (string) "1") {
        dram_nodes[i] = true;
        dram_nodes[i+n] = true;
      }
    }
  } else {
    parallel_for(long i=0; i<n; i++) {
      //std::cout << "n " << i << " ";
      if (dn.Strings[i] == (string) "1") {
      //std::cout << "string " ;
        dram_nodes[i] = true;
      //std::cout << "i " ;
        dram_nodes[i+n] = true;
      }
    }
  }

  std::cout << "done dram nodes\n";

  uintE dram_count = 0;
  uintE other_count = 0;
  // put the indegrees together and the outdegrees together
  // rather than the indegrees and outdegress of a node together
  for (size_t i=0; i<n; i++) {
    if (dram_nodes[i]) {
      partition_offsets[i] = dram_count;
      dram_count += ((i == n-1) ? m : o[i+1])-o[i];
    } else {
      partition_offsets[i] = other_count;
      other_count += ((i == n-1) ? m : o[i+1])-o[i];
    }
  }

  std::cout << "done outgoing partition offsets\n";
  for (size_t i=0; i<n; i++) {
    if (dram_nodes[i]) {
      partition_offsets[i+n] = dram_count;
      dram_count += indo[i];
    } else {
      partition_offsets[i+n] = other_count;
      other_count += indo[i];
    }
  }

  //for (size_t i=0; i<_n; i++)
  //  std::cout <<"dram node pre |" << dram_nodes[i] << std::endl;
  std::cout << "created node partitioner from file" << std::endl;

}

NodePartitioner::NodePartitioner(size_t n, size_t m, uintE* o, size_t capacity) {
  _n = n;
  dram_nodes = newA(bool, n);
  partition_offsets = newA(uintE, n);

  if (capacity >= m) {
    //std::cerr << "capacity is everything" << std::endl;
    memset(dram_nodes, true, n);
    partition_offsets = o;
    return;
  }

  memset(dram_nodes, false, n);

  // add cilk in place sort later
  std::vector<std::pair<long,uintE>> o_vec; // do 2 in place sorts later
  o_vec.resize(n);

  //std::cout << " o[0] " << o[0] << "," << o[1] << std::endl;

  parallel_for(long i=0; i < n-1; i++) {
    o_vec[i] = std::pair<long,uintE>(i, o[i+1]-o[i]);
  }

  o_vec[n-1] = std::pair<long,uintE>(n-1, m - o[n-1]);

  std::sort(o_vec.begin(), o_vec.end(),
      [](const std::pair<long,uintE>& a, const std::pair<long,uintE>& b) {
        return a.second > b.second;
      });

  size_t i = 0;
  size_t occupied = 0;
  //std::cerr << "occupied <cap? " << occupied << "+" << o_vec[i].second << std::endl;
  while((i<n)){
    if ((occupied+o_vec[i].second) <= capacity) {
    dram_nodes[o_vec[i].first] = true;
    occupied += o_vec[i].second;
    }
    i++;
  }
  //while((occupied+o_vec[i].second) <= capacity) {
  //  dram_nodes[o_vec[i].first] = true;
  //  occupied += o_vec[i].second;
  //  //std::cout << " in nvram => " << o_vec[i].first  << "," << o_vec[i].second << std::endl;
  //  i++;
  //}
  //std::cout << "stopped at i = " << i << std::endl;

  uintE dram_count = 0;
  uintE other_count = 0;
  for (size_t i=0; i<n; i++) {
    if (dram_nodes[i]) {
      partition_offsets[i] = dram_count;
      dram_count += ((i == n-1) ? m : o[i+1])-o[i];
    } else {
      partition_offsets[i] = other_count;
      other_count += ((i == n-1) ? m : o[i+1])-o[i];
    }
  }
  std::cout << "Created node partitioner with capacity " << capacity << std::endl;
  //for (size_t i=0; i<5; i++) { std::cerr << "first 5 dram_node | " << dram_nodes[i] << std::endl; }
  //for (size_t i=0; i<5; i++) { std::cerr << "first 5 partition_offsets | " << partition_offsets[i] << std::endl; }
  //for (size_t i=0; i<n; i++) { std::cerr << "dram_node | " << dram_nodes[i] << std::endl; }
  //for (size_t i=0; i<n; i++) { std::cerr << "partition_offsets | " << partition_offsets[i] << std::endl; }
}

NodePartitioner::NodePartitioner(size_t n, size_t m, uintE* o, uintT* indo, size_t capacity) {
  _n = n;
  dram_nodes = newA(bool, 2*n);
  partition_offsets = newA(uintE, 2*n);

  std::cout << "Using asymmetric graph...capacity is " << capacity << " out of " << 2*m << std::endl;

  if (capacity >= 2*m) {
    //std::cerr << "capacity is everything" << std::endl;
    memset(dram_nodes, true, 2*n);
    memcpy(partition_offsets,o,n*sizeof(uintE));
  } else {

    memset(dram_nodes, false, 2*n);

    // add cilk in place sort later
    std::vector<std::pair<long,uintE>> o_vec; // do 2 in place sorts later
    o_vec.resize(2*n);


    parallel_for(long i=0; i < n-1; i++) {
      o_vec[i] = std::pair<long,uintE>(i, o[i+1]-o[i]);
      o_vec[i+n]= std::pair<long,uintE>(i+n, indo[i]);
    }
    o_vec[n-1] = std::pair<long,uintE>(n-1, m - o[n-1]);
    o_vec[n+n-1]= std::pair<long,uintE>(n+n-1, indo[n-1]);

    std::sort(o_vec.begin(), o_vec.end(),
        [](const std::pair<long,uintE>& a, const std::pair<long,uintE>& b) {
          return a.second > b.second;
        });

    size_t i = 0;
    size_t occupied = 0;
    while((i<2*n)){
      if ((occupied+o_vec[i].second) <= capacity) {
      dram_nodes[o_vec[i].first] = true;
      occupied += o_vec[i].second;
      }
      i++;
    }
  }
  uintE dram_count = 0;
  uintE other_count = 0;
  for (size_t i=0; i<n; i++) {
    if (dram_nodes[i]) {
      partition_offsets[i] = dram_count;
      dram_count += ((i == n-1) ? m : o[i+1])-o[i];
    } else {
      partition_offsets[i] = other_count;
      other_count += ((i == n-1) ? m : o[i+1])-o[i];
    }
  }
  for (size_t i=0; i<n; i++) {
    if (dram_nodes[n+i]) {
      partition_offsets[n+i] = dram_count;
      dram_count += indo[i];
    } else {
      partition_offsets[n+i] = other_count;
      other_count += indo[i];
    }
  }

  std::cout << "Created node partitioner with capacity " << capacity << std::endl;
  //for (size_t i=0; i<5; i++) { std::cerr << "first 5 dram_node | " << dram_nodes[i] << std::endl; }
  //for (size_t i=0; i<5; i++) { std::cerr << "first 5 partition_offsets | " << partition_offsets[i] << std::endl; }
  //for (size_t i=0; i<2*n; i++) { std::cerr << "dram_node | " << dram_nodes[i] << std::endl; }
  //for (size_t i=0; i<2*n; i++) { std::cerr << "partition_offsets | " << partition_offsets[i] << std::endl; }
}

NodePartitioner::~NodePartitioner() {
  if (dram_nodes != nullptr) free(dram_nodes);
  if (partition_offsets != nullptr) free(partition_offsets);
}
#endif //NODE_PARTITIONER
