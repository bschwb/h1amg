#ifndef MYAMG_SAMPLE_SORT_HPP_
#define MYAMG_SAMPLE_SORT_HPP_

#include <random>
#include <iterator>

#include <ngstd.hpp>

namespace h1amg
{

// Adapted from http://www.cppsamples.com/common-tasks/choose-random-element.html
int select_randomly(FlatArray<int> index) {
  static std::random_device rd;
  static std::mt19937 engine(rd());
  std::uniform_int_distribution<int> dist(0, index.Size()-1);

  return index[dist(engine)];
}


void SampleSortI(FlatArray<double> data, FlatArray<int> index)
{
  static Timer Tsample_sort("Sample Sort");
  RegionTimer Rsample_sort(Tsample_sort);

  size_t n = index.Size();
  int nr_buckets = ngstd::TaskManager::GetNumThreads();
  int sample_size = nr_buckets-1;
  int over_sample = 10;
  int over_sample_size = over_sample * sample_size;

  static Timer T1("Sample Sort - over sample pivots");
  T1.Start();
  Array<int> over_sampled_ind(over_sample_size);
  for (auto i : Range(over_sample_size)) {
    over_sampled_ind[i] = select_randomly(index);
  }
  T1.Stop();

  static Timer T2("Sample Sort - sort over sampled pivots");
  T2.Start();
  QuickSortI(data, over_sampled_ind);
  T2.Stop();

  static Timer T2b("Sample Sort - select pivots");
  T2b.Start();
  Array<int> sample_ind(sample_size);
  for (auto i : Range(sample_size)) {
    sample_ind[i] = over_sampled_ind[i * over_sample];
  }
  T2b.Stop();

  static Timer T3("Sample Sort - place in buckets");
  T3.Start();

  static Timer T3_1("Sample Sort - determine bucket");
  T3_1.Start();
  Array<int> bucket_of_ind(n);
  ParallelFor(n, [&] (auto i) {
      int start = 0;
      int end = sample_size-1;
      int mid = (start+end)/2;
      while (start <= end) {
        mid = (start+end)/2;
        if (data[sample_ind[mid]] < data[i]) {
          start = mid + 1;
          continue;
        }
        else if (data[sample_ind[mid]] > data[i]) {
          end = mid - 1;
          continue;
        }
        else { break; }
      }
      if (start > end) {
        bucket_of_ind[i] = start;
      }
      else {
        bucket_of_ind[i] = mid;
      }
    });
  T3_1.Stop();

  static Timer T3_2("Sample Sort - inverse index bucket map");
  T3_2.Start();
  ngstd::TableCreator<int> buckets_creator(nr_buckets);
  /*
  for (; !buckets_creator.Done(); buckets_creator++) {
    ParallelFor(nr_buckets, [&] (auto bucket) {
      for (auto i : Range(n)) {
        if (bucket == bucket_of_ind[i]) {
          buckets_creator.Add(bucket, i);
        }
      }
    });
  }
  */
  /*
  for (; !buckets_creator.Done(); buckets_creator++) {
    ParallelFor (n, [&] (auto i)
                 {
                   buckets_creator.Add (bucket_of_ind[i], i);
                 });
  }
  */

  for (; !buckets_creator.Done(); buckets_creator++) {
    ParallelForRange (n, [&] (IntRange r)
                      {
                        ngstd::TableCreator<int> mycreator(nr_buckets);
                        for (; !mycreator.Done(); mycreator++)                      
                          for (auto i : r)
                            mycreator.Add (bucket_of_ind[i], i);

                        auto mytab = mycreator.MoveTable();
                        for (auto i : Range(nr_buckets))
                          buckets_creator.Add (i, mytab[i]);
                      });
  }
  
  auto table = buckets_creator.MoveTable();
  T3_2.Stop();
  T3.Stop();

  static Timer T4("Sample Sort - sort buckets");
  T4.Start();
  ParallelFor(nr_buckets, [&] (auto bucket) {
    QuickSortI(data, table[bucket]);
  });
  T4.Stop();

  static Timer T5("Sample Sort - merge indices");
  T5.Start();
  int start = 0;
  int end = 0;
  for (int bucket = 0; bucket < table.Size(); ++bucket) {
    end += table[bucket].Size();
    index.Range(start, end) = table[bucket];
    start = end;
  }
  T5.Stop();
  /*
  for (size_t i = 0; i < n-1; i++)
    if (data[index[i]] > data[index[i+1]])
      cout << "sort wrong" << endl;
  */
}


}  // myamg

#endif  // MYAMG_SAMPLE_SORT_HPP_
