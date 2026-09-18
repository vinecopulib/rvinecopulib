// Copyright © 2016-2026 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <chrono>
#include <condition_variable>
#include <exception>
#include <functional>
#include <future>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

namespace vinecopulib {

namespace tools_thread {

//! Implementation of the thread pool pattern based on `std::thread`.
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
  // `jobs_`, `stopped_`, `num_busy_` and `error_ptr_` may only be touched
  // while holding `m_tasks_`. A member whose name ends in `_locked` requires
  // the caller to hold it and keeps it held; one taking a
  // `std::unique_lock<std::mutex>&` requires the caller to hold it, but may
  // release it while waiting.
  void start_worker();
  void do_job(std::function<void()>&& job);
  void announce_busy_locked();
  void announce_idle();
  void announce_stop();
  void join_workers();
  void clear_locked();

  bool has_errored_locked() const;
  bool all_jobs_done_locked() const;
  void wait_for_jobs();
  bool wait_for_wake_up_event(std::unique_lock<std::mutex>& lk);
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
{
}

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
  if (workers_.empty()) {
    f(args...); // if there are no workers, do the job in the main thread
    return;
  } else {
    // must hold lock while modifying the shared queue
    std::lock_guard<std::mutex> lk(m_tasks_);
    if (stopped_)
      throw std::runtime_error("cannot push to joined thread pool");
    // bind moves/copies the decayed arguments instead of capturing them by
    // value a second time
    jobs_.emplace(std::bind(std::forward<F>(f), std::forward<Args>(args)...));
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

//! @brief Waits for all jobs to finish, but does not join the threads.
//!
//! @details A job's exception is rethrown here once, and cancels the jobs that
//! have not started. The pool stays usable afterwards.
inline void
ThreadPool::wait()
{
  this->wait_for_jobs();
  this->rethrow_exceptions();
}

//! @brief Waits for all jobs to finish and joins all threads.
//!
//! @details The threads are stopped and joined even when a job threw.
inline void
ThreadPool::join()
{
  this->wait_for_jobs();
  this->announce_stop();
  this->join_workers();
  this->rethrow_exceptions();
}

//! clears the pool from all open jobs.
inline void
ThreadPool::clear()
{
  // must hold lock while modifying job queue
  std::lock_guard<std::mutex> lk(m_tasks_);
  this->clear_locked();
}

//! clears the pool from all open jobs (must be called while holding
//! `m_tasks_`).
inline void
ThreadPool::clear_locked()
{
  std::queue<std::function<void()>>().swap(jobs_);
  cv_tasks_.notify_all();
}

//! spawns a worker thread waiting for jobs to arrive.
inline void
ThreadPool::start_worker()
{
  workers_.emplace_back([this] {
    while (true) {
      // must hold a lock while reading or modifying shared variables
      std::unique_lock<std::mutex> lk(m_tasks_);

      // thread should wait when there is no job
      cv_tasks_.wait(lk, [this] { return stopped_ || !jobs_.empty(); });

      // an empty queue implies the pool was stopped, and nothing can be
      // pushed to a stopped pool; there is no work left to wait for
      if (jobs_.empty())
        return;

      // take job from the queue
      auto job = std::move(jobs_.front());
      jobs_.pop();

      // lock can be released before starting work, but must signal
      // that thread will be busy before (!) to avoid premature breaks
      this->announce_busy_locked();
      lk.unlock();

      this->do_job(std::move(job));
      this->announce_idle();
      std::this_thread::yield();
    }
  });
}

//! executes a job safely and let's pool know when it's busy.
//! @param job Job to be executed.
inline void
ThreadPool::do_job(std::function<void()>&& job)
{
  try {
    job();
  } catch (...) {
    {
      std::lock_guard<std::mutex> lk(m_tasks_);
      // the first failure is the one that cancels the queue
      if (!this->has_errored_locked())
        error_ptr_ = std::current_exception();
    }
    cv_busy_.notify_one();
  }
}

//! signals that a worker is busy (must be called while holding `m_tasks_`).
inline void
ThreadPool::announce_busy_locked()
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
  if (!workers_.empty()) {
    for (auto& worker : workers_) {
      if (worker.joinable())
        worker.join();
    }
  }
}

//! checks if an error occurred (must be called while holding `m_tasks_`).
inline bool
ThreadPool::has_errored_locked() const
{
  return static_cast<bool>(error_ptr_);
}

//! check whether all jobs are done (must be called while holding `m_tasks_`).
inline bool
ThreadPool::all_jobs_done_locked() const
{
  return (num_busy_ == 0) && jobs_.empty();
}

//! @brief Waits until no job is queued or running.
inline void
ThreadPool::wait_for_jobs()
{
  // must hold the lock while reading the shared state; the wake up event
  // releases it while waiting
  std::unique_lock<std::mutex> lk(m_tasks_);
  while (true) {
    if (this->wait_for_wake_up_event(lk)) {
      if (this->all_jobs_done_locked())
        return;
      // an error makes the jobs that have not started pointless; the ones
      // already running still have to finish
      this->clear_locked();
    }
  }
}

//! checks whether `wait()` needs to wake up, i.e., all jobs are done or an
//! error makes the jobs that have not started pointless.
//! @param lk A lock on `m_tasks_` held by the caller; released while waiting.
inline bool
ThreadPool::wait_for_wake_up_event(std::unique_lock<std::mutex>& lk)
{
  static auto timeout = std::chrono::milliseconds(250);
  auto wake_up_event_occurred = [this] {
    return this->all_jobs_done_locked() ||
           (this->has_errored_locked() && !jobs_.empty());
  };
  // the timeout bounds the wait: `cv_busy_` is notified to a single waiter,
  // and pushing a job does not notify it at all
  cv_busy_.wait_for(lk, timeout, wake_up_event_occurred);
  return wake_up_event_occurred();
}

//! @brief Rethrows the exception stored by a failing job, and consumes it.
inline void
ThreadPool::rethrow_exceptions()
{
  std::exception_ptr error_ptr;
  {
    // must hold the lock while reading the stored exception
    std::lock_guard<std::mutex> lk(m_tasks_);
    std::swap(error_ptr, error_ptr_);
  }
  if (error_ptr)
    std::rethrow_exception(error_ptr);
}

}
}
