#ifndef DYNPQ_H
#define DYNPQ_H 1
#include <vector>
#include<algorithm>
using namespace std;
// compares the associated values of passed iterators
template<class T>
struct IterGreater {
    bool operator()( T x,  T y) const { return *y < *x;}
};
template <class key_type>
class  dynamic_priority_queue {
  public:
    // public type definitions
    typedef typename std::vector<key_type>::size_type size_type;
    typedef typename std::vector<key_type>::difference_type index_type;

    // constructor
    dynamic_priority_queue(std::vector<key_type>& v);

    // change a value at position 'at'
    void changeKeyAt(index_type at, key_type k);

    // index of the smallest element (= highest priority)
    index_type topIndex() const { return c.front() - first; }

    // value of the smallest element (= highest priority)
    const key_type& topKey() const { return *c.front(); }

    void pop();       // remove smallest element from the heap

    bool empty() const { return csize == 0;}
    size_type size() const { return csize;}

  private:
    //checkedvector<index_type> Indices;      // auxiliary vector
    vector<index_type> Indices;      // auxiliary vector
    typedef typename std::vector<key_type>::iterator randomAccessIterator;
    //checkedvector<randomAccessIterator> c;  // heap of iterators
    vector<randomAccessIterator> c;  // heap of iterators
    randomAccessIterator first;             // beginning of the external vector
    IterGreater<randomAccessIterator> comp; // comparison object
    index_type csize;                       // current heap size

    // heap update (see below)
    void goUp(index_type);
    void goDown(index_type);
};

template <class key_type>
dynamic_priority_queue<key_type>::dynamic_priority_queue(
     std::vector<key_type>& v)
  : Indices(v.size()), c(v.size()), first(v.begin()),
  csize(v.size()) {
    // store iterators and generate heap
    index_type i;
	for(i = 0; i< csize; ++i) 
       c[i] = v.begin()+i;
    make_heap(c.begin(), c.end(), comp);        // STL

    // construct index array
    for(i = 0; i< csize; ++i)
       Indices[c[i] - first] = static_cast<int>(i);
}

template <class key_type>
void dynamic_priority_queue<key_type>::changeKeyAt(
                                index_type at, key_type k) {
   index_type idx = Indices[at];
//   assert(idx < csize);   // value still present in the queue?

   if(*c[idx] != k)       // in case of equality, do nothing
      if(k > *c[idx]) {
           *c[idx] = k;   // enter heavier value
           goDown(idx);   // reorganize heap
      }
      else {
           *c[idx] = k;   // enter lighter value
           goUp(idx);     // reorganize heap
      }
}

template <class key_type>
void dynamic_priority_queue<key_type>::goUp(index_type idx) {
    index_type Predecessor = (idx-1)/2;
    randomAccessIterator temp = c[idx];

    while(Predecessor != idx && comp(c[Predecessor], temp)) {
       c[idx] = c[Predecessor];
       Indices[c[idx]-first] = static_cast<int>(idx);
       idx = Predecessor;
       Predecessor = (idx-1)/2;
    }
    
    c[idx] = temp;
    Indices[c[idx]-first] = static_cast<int>(idx);
}

template <class key_type>
void dynamic_priority_queue<key_type>::goDown(index_type idx) {
     index_type Successor = (idx+1)*2-1;

     if(Successor < csize-1
        && comp(c[Successor], c[Successor+1]))
           ++Successor;
     randomAccessIterator temp = c[idx];

     while(Successor < csize && comp(temp, c[Successor])) {
        c[idx] = c[Successor];
        Indices[c[idx]-first] = static_cast<int>(idx);
        idx = Successor;
        Successor = (idx+1)*2-1;

        if(Successor < csize-1
           && comp(c[Successor], c[Successor+1]))
              ++Successor;
     }
     c[idx] = temp;
     Indices[c[idx]-first] = static_cast<int>(idx);
}

/* The method pop() removes the topmost element from the heap. This is
   done by moving the last element to the top and blocking the freed
   position with --csize. Subsequently, the element sinks down to its
   proper position. */

template <class key_type>
void dynamic_priority_queue<key_type>::pop() {
   // overwrite iterator at the top with the 
   // address of the last element
   c[0] = c[--csize];

   // enter the new address 0 at the position belonging 
   // to this element in the auxiliary array
   Indices[c[0]-first] = 0;

   // let the element at the top sink to the correct
   // position corresponding to its weight
   goDown(0);
}

#endif





