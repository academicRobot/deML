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
    : q(), m(), c(), max_size(max_size)
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
      c.wait(lock);
    }
    q.push(t);
    c.notify_one();
  }

  // Get the "front"-element.
  // If the queue is empty, wait till a element is avaiable.
  T dequeue(void)
  {
    std::unique_lock<std::mutex> lock(m);
    while(q.empty())
    {
      // release lock as long as the wait and reaquire it afterwards.
      c.wait(lock);
    }
    T val = q.front();
    q.pop();
    return val;
  }

private:
  std::queue<T> q;
  mutable std::mutex m;
  std::condition_variable c;
  size_t max_size;
};

#endif
