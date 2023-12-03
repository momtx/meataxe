////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Parallel execution (threads) support
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>

/// @private
struct PexGroup {
   struct PexGroup* next;
   struct PexGroup** prev;
   uint32_t groupId;
   void* userData;
   void (*finalize)(void* userData);
   #if defined(MTX_DEFAULT_THREADS)
   pthread_mutex_t mutex;
   #endif
   size_t nPending;             // number of tasks waiting or being executed
};

#if defined(MTX_DEFAULT_THREADS)
static pthread_mutex_t groupsMutex = PTHREAD_MUTEX_INITIALIZER;
#endif
static PexGroup_t* groupsHead = NULL;
static PexGroup_t** groupsTail = &groupsHead;

/// @private
typedef struct Task {
   struct Task* next;
   struct PexGroup* group;
   void (*f)(void* userData);
   void (*fr)(void* userData, size_t begin, size_t end);
   void* userData;
   size_t begin;
   size_t end;
} Task_t;


/// @private
struct ThreadInfo {
   #if defined(MTX_DEFAULT_THREADS)
   pthread_t id;
   #endif
   int threadNumber;    // 1...mThreads (0 is reserved for the main thread)
   char name[100];
   char logPrefix[105];
   struct ErrorContextStack contextStack;
};

static int tidWidth = 0;

#if defined(MTX_DEFAULT_THREADS)

static int threadPoolSize = 0;        // Fixed (from -j)
static pthread_key_t tidKey;
static pthread_t* threadId;           // Worker threads only.

#else

static const int threadPoolSize = 0;

#endif

#if defined(MTX_DEFAULT_THREADS)
static pthread_once_t createTidKeyOnce = PTHREAD_ONCE_INIT;
#endif

// Data protected by tqMutex
#if defined(MTX_DEFAULT_THREADS)
static pthread_mutex_t tqMutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t tqWakeup = PTHREAD_COND_INITIALIZER;   // task added or shutdown started
static pthread_cond_t tqIdle = PTHREAD_COND_INITIALIZER;     // task finished
static int nThreads = 0;           // Number of created threads
static int nBusyThreads = 0;       // Number of threads executing a task
static int tqShutdown = 0;
static struct Task* tqHead = NULL;
static struct Task** tqTail = &tqHead;
#endif
static int groupId = 0;

////////////////////////////////////////////////////////////////////////////////////////////////////

#if defined(MTX_DEFAULT_THREADS)
static void deleteGroup(PexGroup_t* group)
{
   MTX_ASSERT(group->prev != NULL);
   pthread_mutex_lock(&groupsMutex);
   if ((*group->prev = group->next) == NULL)
      groupsTail = group->prev;
   else
      group->next->prev = group->prev;
   pthread_mutex_unlock(&groupsMutex);
   pthread_mutex_destroy(&group->mutex);
   memset(group, 0, sizeof(PexGroup_t));
   sysFree(group);
}
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////

static void setThreadName(struct ThreadInfo* ti, const char* name)
{
   snprintf(ti->name, sizeof(ti->name), "%s", name);
   if (*name == 0) {
      snprintf(ti->logPrefix, sizeof(ti->logPrefix), "[%*d] ", tidWidth, ti->threadNumber);
   }
   else {
      snprintf(ti->logPrefix, sizeof(ti->logPrefix), "[%*d:%s] ", tidWidth, ti->threadNumber, name);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Creates a task group.
/// Task groups are used to execute a "finalizer" function after all tasks in the group are done.
/// Using a task group always requires the following steps:
/// 1. Create the task group.
/// 2. Add arbitrary tasks to the group by passing the group as first argument to @ref pexExecute.
/// 3. Call @ref pexFinally to define the finalizer.
///
/// It is an error to create a task group without calling @ref pexFinally for this
/// group. The error will be detected by @ref pexShutdown.

PexGroup_t* pexCreateGroup()
{
   PexGroup_t* group = ALLOC(PexGroup_t);
   memset(group, 0, sizeof(PexGroup_t));
   #if defined(MTX_DEFAULT_THREADS)
   pthread_mutex_init(&group->mutex, NULL);
   pthread_mutex_lock(&groupsMutex);
   #endif
   group->groupId = ++groupId;
   MTX_ASSERT(groupsTail != NULL);
   MTX_ASSERT(*groupsTail == NULL);
   group->prev = groupsTail;
   *groupsTail = group;
   groupsTail = &group->next;
   #if defined(MTX_DEFAULT_THREADS)
   pthread_mutex_unlock(&groupsMutex);
   #endif
   return group;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

#if defined(MTX_DEFAULT_THREADS)
static void removeTaskFromGroup(PexGroup_t* group)
{
   pthread_mutex_lock(&group->mutex);
   MTX_ASSERT(group->nPending > 0);
   --group->nPending;
   const int finalize = group->nPending == 0 && group->finalize != NULL;
   MTX_LOG2("%s: nPending=%lu, finalize=%p",
         __func__, (unsigned long) group->nPending, group->finalize);
   pthread_mutex_unlock(&group->mutex);
   if (finalize) {
      group->finalize(group->userData);
      deleteGroup(group);
   }
}
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////

#if defined(MTX_DEFAULT_THREADS)
static void addTaskToGroup(PexGroup_t* group)
{
   pthread_mutex_lock(&group->mutex);
   MTX_ASSERT(group->finalize == NULL);
   ++group->nPending;
   pthread_mutex_unlock(&group->mutex);
}
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Sets the task group finalizer.
/// If PEX is not initialized, pexFinally() calls the finalizer function, @a f, with the argument
/// @a userData and ignores @a group.
///
/// Otherwise the finalizer call is scheduled for execution after all tasks in the given group
/// are completed. If this is already the case, the finalizer is executed immediately in the 
/// calling thread.
///
/// pexFinally() must be called exactly once for any group. Calling pexFinally() a second time
/// fails immediately. Task groups without finalizers are detected in @ref pexShutdown and cause
/// an error, too.
///
/// After return from pexFinally() the group pointer must be treated as invalid and not be passed
/// to any other PEX function. In particular, no more tasks can be added to the group.
///
/// pexFinally() may be called from any thread. In particular, pexFinally() may be called from a
/// task in the same group.
///

void pexFinally(PexGroup_t *group, void (*f)(void* userData), void* userData)
{
   MTX_ASSERT(group != NULL);
   MTX_ASSERT(f != NULL);
   if (threadPoolSize == 0) {
      // Execute immediately.
      f(userData);
      return;
   }

#if defined(MTX_DEFAULT_THREADS)
   pthread_mutex_lock(&group->mutex);
   MTX_ASSERT(group->finalize == NULL);
   MTX_LOG2("%s: grp=%p nPending=%lu", __func__, group, (unsigned long)group->nPending);
   if (group->nPending != 0) {
      // schedule for execution after last task
      group->finalize = f;
      group->userData = userData;
      pthread_mutex_unlock(&group->mutex);
   } else {
      // finalize immediately
      pthread_mutex_unlock(&group->mutex);
      f(userData);
      deleteGroup(group);
   }
#endif
}

////////////////////////////////////////////////////////////////////////////////////////////////////

#if defined(MTX_DEFAULT_THREADS)

static void executeTask(Task_t* task)
{
   const int context = mtxBegin(MTX_HERE, "executing task");
   const unsigned long groupId = (task->group != NULL) ? task->group->groupId : 0;
   if (task->fr != NULL) {
      MTX_LOG2("begin task %lx.%p.%p [%lu,%lu)",
         groupId, task->userData, task->fr, (unsigned long) task->begin, (unsigned long) task->end);
      task->fr(task->userData, task->begin, task->end);
   }
   else {
      MTX_LOG2("begin task %lx.%p.%p", groupId, task->userData, task->f);
      task->f(task->userData);
   }
   MTX_LOG2("end task %p.%p", task->userData, task->fr ? (void*)task->fr : (void*)task->f);
   if (task->group) {
      removeTaskFromGroup(task->group);
   }
   mtxEnd(context);
}

#endif

////////////////////////////////////////////////////////////////////////////////////////////////////

#if defined(MTX_DEFAULT_THREADS)

static void createTidKey()
{
    pthread_key_create(&tidKey, NULL);
}

#endif

////////////////////////////////////////////////////////////////////////////////////////////////////

static struct ThreadInfo* getThreadInfo()
{
   #if defined(MTX_DEFAULT_THREADS)

   pthread_once(&createTidKeyOnce, createTidKey);
   struct ThreadInfo* ti = pthread_getspecific(tidKey);
   if (ti == NULL) {
      ti = ALLOC(struct ThreadInfo);
      pthread_setspecific(tidKey, ti);
   }
   // TODO: LINK (using cstackMutex)
   return ti;

   #else

   static struct ThreadInfo mainThreadInfo = {0};
   // TODO: LINK (using cstackMutex)
   return &mainThreadInfo;

   #endif
}

#if defined(MTX_DEFAULT_THREADS)

static void destroyThreadInfo(struct ThreadInfo* ti)
{
   // TODO: UNLINK (using cstackMutex)
   sysFree(ti);
   pthread_setspecific(tidKey, NULL);
}

#endif

////////////////////////////////////////////////////////////////////////////////////////////////////

//void pexForEach(void (*f)(int threadNumber, struct ContextStack* cs), int skipThisThread)
//{
//   // TODO: synchronize with csMutex
//   for (struct ErrorContextStack** pcs = &contextStacks; *pcs; pcs = &(*pcs)->next_) {
//      if (skipThisThread && pthread_equal(getThreadInfo->id, pthread_self())) {
//         continue;
//      }
//      f(cs
//   }
//}

////////////////////////////////////////////////////////////////////////////////////////////////////

#if defined(MTX_DEFAULT_THREADS)

static void threadCleanup()
{
   struct ThreadInfo *ti = getThreadInfo();
   destroyThreadInfo(ti);
}

#endif

////////////////////////////////////////////////////////////////////////////////////////////////////

#if defined(MTX_DEFAULT_THREADS)

static void threadInit(void* arg)
{
   int *threadNumber = (int*) arg;
   struct ThreadInfo *ti = getThreadInfo();
   ti->threadNumber = *threadNumber;
   sysFree(threadNumber);
   setThreadName(ti, "");
}

#endif 

////////////////////////////////////////////////////////////////////////////////////////////////////

#if defined(MTX_DEFAULT_THREADS)

static void* threadMain(void* arg)
{
   pthread_mutex_lock(&tqMutex);
   threadInit(arg);
   MTX_LOG2("worker thread ready");
   while (1) {
       if (tqShutdown)
	   break;
       else if (tqHead != NULL) {
	   struct Task* task = tqHead;
	   if ((tqHead = tqHead->next) == NULL) {
	       tqTail = &tqHead;
	   }
           ++nBusyThreads;
           pthread_mutex_unlock(&tqMutex);
	   task->next = NULL;
           setThreadName(getThreadInfo(), "");
	   executeTask(task);
           pthread_mutex_lock(&tqMutex);
           --nBusyThreads;
           pthread_cond_broadcast(&tqIdle);
	   sysFree(task);
       } else {
          pthread_cond_wait(&tqWakeup, &tqMutex);
       }
   }
   MTX_LOG2("worker thread exiting");
   threadCleanup();
   pthread_mutex_unlock(&tqMutex);
   return NULL;
}

#endif

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns a unique number in the range [0,poolSize] which identifies thge calling thread.
/// Thread number 0 is the program's main thread, numbers greater than 0 are used for member of the
/// thread pool.

int pexThreadNumber()
{
   return getThreadInfo()->threadNumber;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns a short text to be used in log messages to identify the logging thred. By default, the
/// prefix is "[N]", where N ist the number returned by @ref pexThreadNumber. If a thread name is
/// set (see @ref pexSetThreadName), the prefix becomes "[N:NAME]".

const char* pexLogPrefix()
{
   return getThreadInfo()->logPrefix;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Returns the thread name. The thread name is always empty when a thread starts working on a task.
/// A new thread name, which need not be unique, can be set at any time with @ref pexSetThreadName.
const char* pexThreadName()
{
#if defined(MTX_DEFAULT_THREADS)
   struct ThreadInfo* ti = (struct ThreadInfo*) pthread_getspecific(tidKey);
   return ti ? ti->name : "";
#else
   return "";
#endif
}

////////////////////////////////////////////////////////////////////////////////////////////////////

struct ErrorContextStack* pexContextStack()
{
   return &getThreadInfo()->contextStack;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Sets a short name for the thread, which will be used in log messages.
/// The thread name is always reset to "" when a thread starts working on a task.

void pexSetThreadName(const char* name, ...)
{
   struct ThreadInfo* ti = getThreadInfo();
   va_list args;
   va_start(args, name);
   char formattedName[100];
   vsnprintf(formattedName, sizeof(formattedName), name, args);
   va_end(args);
   setThreadName(ti, formattedName);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Waits until all pending tasks are finished.

// TODO: pexWait(PexGroup* group) ?

void pexWait()
{
#if defined(MTX_DEFAULT_THREADS)
   pthread_mutex_lock(&tqMutex);
   while (1) {
      if (tqHead == NULL && nBusyThreads == 0)
         break;
      MTX_LOG2("Waiting for pending tasks");
      pthread_cond_wait(&tqIdle, &tqMutex);
   }
   pthread_mutex_unlock(&tqMutex);
#endif
}


////////////////////////////////////////////////////////////////////////////////////////////////////

/// Stops worker threads and releases internal resources.
/// This function fails if it is called while tasks are pending or being executed.
/// Calling pexShutdown() multiple times or without a prior pexInit() call is allowed and has no
/// effect.

void pexShutdown()
{
   // Allow multiple calls
   if (threadPoolSize == 0)
      return;
      
#if defined(MTX_DEFAULT_THREADS)
   MTX_LOG2("PEX shutting down");

   pthread_mutex_lock(&tqMutex);
   MTX_ASSERT(tqHead == NULL);
   MTX_ASSERT(tqTail == &tqHead);
   MTX_ASSERT(nBusyThreads == 0);
   tqShutdown = 1;
   pthread_cond_broadcast(&tqWakeup);
   pthread_mutex_unlock(&tqMutex);
   for (int i = 0; i < nThreads; ++i) {
      pthread_join(threadId[i], NULL);
   }
   pthread_mutex_unlock(&tqMutex);
   // TODO: REALLOC(threadInfo, struct ThreadInfo, 1);
   sysFree(threadId);
   threadId = NULL;
   threadPoolSize = 0;
   nThreads = 0;

   pthread_mutex_lock(&groupsMutex);
   MTX_ASSERT(groupsHead == NULL);
   pthread_mutex_unlock(&groupsMutex);
#endif
   groupId = 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Initializes the parallel execution (PEX) library.
/// The argument determines the number of worker threads, i.e., the maximum number of tasks that
/// can be executed concurrently, not counting the program's main thread.
/// 
/// This function will fail if PEX is already initialized. To change the thread pool size, call
/// @ref pexShutdown first, then reinitialize the library with the new number of threads.

void pexInit(int poolSize)
{
   MTX_ASSERT(poolSize > 0);
   if (threadId != NULL)
      mtxAbort(MTX_HERE, "Multiple calls of pexInit()");
   #if defined(MTX_DEFAULT_THREADS)
   MTX_LOG2("PEX initializing, poolSize=%d", poolSize);
   threadPoolSize = poolSize;
   tqShutdown = 0;
   nThreads = 0;
   nBusyThreads = 0;
   threadId = NALLOC(pthread_t, threadPoolSize);
   #endif
   tidWidth = threadPoolSize < 10 ? 1 : (threadPoolSize < 100 ? 2 : (threadPoolSize < 1000 ? 3 : 4));
   setThreadName(getThreadInfo(), "");
}


////////////////////////////////////////////////////////////////////////////////////////////////////

void pexSleep(unsigned timeInMs)
{
#if defined(MTX_DEFAULT_THREADS)
    struct timespec expTime;
    clock_gettime(CLOCK_REALTIME, &expTime);
    expTime.tv_sec += timeInMs / 1000;
    if ((expTime.tv_nsec += (timeInMs % 1000L) * 1000000) > 1000000000L) {
        expTime.tv_nsec -= 1000000000;
        ++expTime.tv_sec;
    }

    static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t cond = PTHREAD_COND_INITIALIZER;
    pthread_mutex_lock(&mutex);
    pthread_cond_timedwait(&cond, &mutex, &expTime);
    pthread_mutex_unlock(&mutex);
#endif
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void createTask(
   PexGroup_t* group,
   void (* f)(void* userData),
   void (* fr)(void* userData, size_t begin, size_t end),
   void* userData,
   size_t begin,
   size_t end)
{
   if (threadPoolSize == 0) {
      // Execute immediately.
      f ? f(userData) : fr(userData, begin, end);
      return;
   }

#if defined(MTX_DEFAULT_THREADS)
   if (group != NULL) {
      addTaskToGroup(group);
   }
   Task_t* task = ALLOC(Task_t);
   task->group = group;
   task->f = f;
   task->fr = fr;
   task->userData = userData;
   task->begin = begin;
   task->end = end;
   const unsigned long gid = (group != NULL) ? group->groupId : 0;
   if (fr == NULL) {
      MTX_LOG2("create task %lu.%p.%p", gid, userData, fr ? (void*)fr : (void*)f);
   }
   else {
      MTX_LOG2("create task %lu.%p.%p [%lu,%lu)",
         gid, userData, fr ? (void*)fr : (void*)f,
         (unsigned long) begin, (unsigned long) end);
   }

   pthread_mutex_lock(&tqMutex);
   *tqTail = task;
   task->next = NULL;
   tqTail = &task->next;

   // Create additional worker threads if necessary and possible.
   if (nBusyThreads == nThreads && nThreads < threadPoolSize) {
      int* threadNumber = ALLOC(int); // memory is released in worker thread
      *threadNumber = nThreads + 1;
      int threadCreateResult = pthread_create(threadId + nThreads, NULL, threadMain, threadNumber);
      MTX_ASSERT(threadCreateResult == 0);
      ++nThreads;
   }

   pthread_cond_signal(&tqWakeup);  // could this deadlock? broadcast required???
   pthread_mutex_unlock(&tqMutex);
#endif
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Schedules a task for execution.
/// If PEX is not initialized, this function calls the function @a f with the argument @a userData
/// and ignores @a group.
/// Otherwise the function call is scheduled to be executed by a worker thread at some point in the
/// future.
///
/// Task functions must not make any assumption on the order of execution. For example, if task B
/// was created after task A, B may already be finished before A is started. If task B is created
/// in task A, B may be finished before A is finished.
///
/// If @a group is not NULL, it must be a pointer to a task group created by @ref pexCreateGroup.
/// The task becomes a member of this group, i.e., the group's finalizer function
/// will not be executed before the task is finished. See @ref pexFinally.
/// pexExecute() fails if @ref pexFinally has already been called on the given group.

void pexExecute(PexGroup_t* group, void (*f)(void *userData), void* userData)
{
   MTX_ASSERT(f != NULL);
   createTask(group, f, NULL, userData, 0, 0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Schedules a task for execution.
/// This function works like @ref pexExecute but the task function receives three arguments,
/// @a userData, @a begin and @a end.

void pexExecuteRange(
      PexGroup_t* group,
      void (*f)(void *userData, size_t begin, size_t end),
      void* userData,
      size_t begin,
      size_t end)
{
   MTX_ASSERT(f != NULL);
   createTask(group, NULL, f, userData, begin, end);
}


// vim:sw=3:et:cin
