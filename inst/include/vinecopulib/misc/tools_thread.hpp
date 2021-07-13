// Copyright © 2016-2021 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <atomic>
#include <condition_variable>
#include <future>
#include <queue>
#include <thread>
#include <vector>

namespace vinecopulib {

namespace tools_thread {

//! Implemenation of the thread pool pattern based on `std::sthread`.
class ThreadPool
{
public:
  ThreadPool(ThreadPool&&) = delete;
  ThreadPool(const ThreadPool&) = delete;
  ThreadPool();
  explicit ThreadPool(size_t nThreads);

  ~ThreadPool() noexcept;

  ThreadPool& operator=(const ThreadPool&) = delete;
  ThreadPool& operator=(ThreadPool&& other) = delete;

  template<class F, class... Args>
  void push(F&& f, Args&&... args);

  template<class F, class I>
  void map(F&& f, I&& items);

  void wait();
  void join();
  void clear();

private:
  void start_worker();
  void do_job(std::function<void()>&& job);
  void announce_busy();
  void announce_idle();
  void announce_stop();
  void join_workers();

  bool has_errored();
  bool all_jobs_done();
  bool wait_for_wake_up_event();
  void rethrow_exceptions();

  std::vector<std::thread> workers_;       // worker threads in the pool
  std::queue<std::function<void()>> jobs_; // the task que

  // variables for synchronization between workers
  std::mutex m_tasks_;
  std::condition_variable cv_tasks_;
  std::condition_variable cv_busy_;
  size_t num_busy_{ 0 };
  bool stopped_{ false };
  std::exception_ptr error_ptr_;
};

//! constructs a thread pool with as many workers as there are cores.
inline ThreadPool::ThreadPool()
  : ThreadPool(std::thread::hardware_concurrency())
{}

//! constructs a thread pool with `nThreads` threads.
//! @param nWorkers Number of worker threads to create; if `nThreads = 0`, all
//!    work pushed to the pool will be done in the main thread.
inline ThreadPool::ThreadPool(size_t nWorkers)
{
  for (size_t w = 0; w < nWorkers; ++w)
    this->start_worker();
}

//! destructor joins all threads if possible.
inline ThreadPool::~ThreadPool() noexcept
{
  // destructors should never throw
  try {
    this->announce_stop();
    this->join_workers();
  } catch (...) {
  }
}

//! pushes jobs to the thread pool.
//! @param f A function taking an arbitrary number of arguments.
//! @param args A comma-seperated list of the other arguments that shall
//!   be passed to `f`.
//!
//! The function returns void; if a job returns a result, use `pushReturn()`.
template<class F, class... Args>
void
ThreadPool::push(F&& f, Args&&... args)
{
  if (workers_.size() == 0) {
    f(args...); // if there are no workers, do the job in the main thread
    return;
  } else {
    // must hold lock while modifying the shared queue
    std::lock_guard<std::mutex> lk(m_tasks_);
    if (stopped_)
      throw std::runtime_error("cannot push to joined thread pool");
    jobs_.emplace([f, args...] { f(args...); });
  }
  // signal a waiting worker that there's a new job
  cv_tasks_.notify_one();
}

//! maps a function on a list of items, possibly running tasks in parallel.
//! @param f Function to be mapped.
//! @param items An objects containing the items on which `f` shall be
//!   mapped; must allow for `auto` loops (i.e., `std::begin(I)`/
//!  `std::end(I)` must be defined).
template<class F, class I>
void
ThreadPool::map(F&& f, I&& items)
{
  for (auto&& item : items)
    this->push(f, item);
}

//! waits for all jobs to finish, but does not join the threads.
inline void
ThreadPool::wait()
{
  while (true) {
    if (wait_for_wake_up_event()) {
      if (!this->all_jobs_done()) {
        this->clear(); // cancel all remaining jobs
        continue;      // wait for currently running jobs
      }
      if (this->has_errored() | this->all_jobs_done())
        break;
    }
    std::this_thread::yield();
  }

  this->rethrow_exceptions();
}

//! waits for all jobs to finish and joins all threads.
inline void
ThreadPool::join()
{
  this->wait();
  this->announce_stop();
  this->join_workers();
}

//! clears the pool from all open jobs.
inline void
ThreadPool::clear()
{
  // must hold lock while modifying job queue
  std::lock_guard<std::mutex> lk(m_tasks_);
  std::queue<std::function<void()>>().swap(jobs_);
  cv_tasks_.notify_all();
}

//! spawns a worker thread waiting for jobs to arrive.
inline void
ThreadPool::start_worker()
{
  workers_.emplace_back([this] {
    std::function<void()> job;
    // observe thread pool; only stop after all jobs are done
    while ((!stopped_) | (!jobs_.empty())) {
      // must hold a lock while modifying shared variables
      std::unique_lock<std::mutex> lk(m_tasks_);

      // thread should wait when there is no job
      cv_tasks_.wait(lk, [this] { return stopped_ || !jobs_.empty(); });

      // queue can be empty if thread pool is stopped
      if (jobs_.empty())
        continue;

      // take job from the queue
      job = std::move(jobs_.front());
      jobs_.pop();

      // lock can be released before starting work, but must signal
      // that thread will be busy before (!) to avoid premature breaks
      this->announce_busy();
      lk.unlock();

      this->do_job(std::move(job));
      this->announce_idle();
      std::this_thread::yield();
    }
  });
}

//! executes a job safely and let's pool know when it's busy.
//! @param job Job to be exectued.
inline void
ThreadPool::do_job(std::function<void()>&& job)
{
  try {
    job();
  } catch (...) {
    std::lock_guard<std::mutex> lk(m_tasks_);
    error_ptr_ = std::current_exception();
  }
}

//! signals that a worker is busy (must be called why locking m_tasks_).
inline void
ThreadPool::announce_busy()
{
  ++num_busy_;
  cv_busy_.notify_one();
}

//! signals that a worker is idle.
inline void
ThreadPool::announce_idle()
{
  {
    std::lock_guard<std::mutex> lk(m_tasks_);
    --num_busy_;
  }
  cv_busy_.notify_one();
}

//! signals threads that no more new work is coming.
inline void
ThreadPool::announce_stop()
{
  {
    std::unique_lock<std::mutex> lk(m_tasks_);
    stopped_ = true;
  }
  cv_tasks_.notify_all();
}

//! joins worker threads if possible.
inline void
ThreadPool::join_workers()
{
  if (workers_.size() > 0) {
    for (auto& worker : workers_) {
      if (worker.joinable())
        worker.join();
    }
  }
}

//! checks if an error occured.
inline bool
ThreadPool::has_errored()
{
  return static_cast<bool>(error_ptr_);
}

//! check whether all jobs are done
inline bool
ThreadPool::all_jobs_done()
{
  return (num_busy_ == 0) && jobs_.empty();
}

//! checks whether `wait()` needs to wake up
inline bool
ThreadPool::wait_for_wake_up_event()
{
  static auto timeout = std::chrono::milliseconds(250);
  auto wake_up_event_occured = [this] {
    return this->all_jobs_done() | this->has_errored();
  };
  std::unique_lock<std::mutex> lk(m_tasks_);
  cv_busy_.wait_for(lk, timeout, wake_up_event_occured);
  return wake_up_event_occured();
}

//! rethrows exceptions (exceptions from workers are caught and stored; the
//! wait loop only checks, but does not throw exceptions)
inline void
ThreadPool::rethrow_exceptions()
{
  if (this->has_errored())
    std::rethrow_exception(error_ptr_);
}

}
}
