#include "index_utility.cpp"

using namespace seqan;


template <typename TDna, unsigned span>
bool createHIndex(StringSet<String<TDna> > & seq, HIndex<span> & index, unsigned & threads, bool efficient)
{
  //  if (threads > 1)
  //  {
        return _createQGramIndexDirSA_parallel(seq, index.xstr, index.ysa, index.shape, index.emptyDir, threads, efficient);
  //  }
  //  else 
  //  {
   //     return _createQGramIndexDirSA(seq, index.xstr, index.ysa, index.shape, index.emptyDir);
   // }
        
}