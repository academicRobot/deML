#ifndef SAFE_QUEUE
#define SAFE_QUEUE

#include <queue>
#include <mutex>
#include <condition_variable>

// A threadsafe-queue. https://stackoverflow.com/a/16075550
template <class T>
class SafeQueue
{
public:
  SafeQueue(size_t max_size)
    : q(), m(), empty_cv(), full_cv(), max_size(max_size)
  {}

  ~SafeQueue(void)
  {}

  // Add an element to the queue.
  void enqueue(T t)
  {
    std::unique_lock<std::mutex> lock(m);
    while(q.size()==max_size)
    {
      // release lock as long as the wait and reaquire it afterwards.
      full_cv.wait(lock);
    }
    q.push(t);
    empty_cv.notify_one();
  }

  // Get the "front"-element.
  // If the queue is empty, wait till a element is avaiable.
  T dequeue(void)
  {
    std::unique_lock<std::mutex> lock(m);
    while(q.empty())
    {
      // release lock as long as the wait and reaquire it afterwards.
      empty_cv.wait(lock);
    }
    T val = q.front();
    q.pop();
    full_cv.notify_one();
    return val;
  }

private:
  std::queue<T> q;
  mutable std::mutex m;
  std::condition_variable empty_cv, full_cv;
  size_t max_size;
};

#endif
