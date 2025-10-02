#define _GNU_SOURCE
#include <dlfcn.h>
#include <malloc.h>  // malloc_usable_size (glibc)
#include <signal.h>
#include <stdatomic.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

static void *(*real_malloc)(size_t) = NULL;
static void (*real_free)(void *) = NULL;
static void *(*real_realloc)(void *, size_t) = NULL;
static void *(*real_calloc)(size_t, size_t) = NULL;
static int (*real_posix_memalign)(void **, size_t, size_t) = NULL;
static void *(*real_aligned_alloc)(size_t, size_t) = NULL;

static atomic_uint_least64_t cur_bytes = 0;
static atomic_uint_least64_t peak_bytes = 0;

static void resolve_syms(void) {
  real_malloc = dlsym(RTLD_NEXT, "malloc");
  real_free = dlsym(RTLD_NEXT, "free");
  real_realloc = dlsym(RTLD_NEXT, "realloc");
  real_calloc = dlsym(RTLD_NEXT, "calloc");
  real_posix_memalign = dlsym(RTLD_NEXT, "posix_memalign");
  real_aligned_alloc =
      dlsym(RTLD_NEXT, "aligned_alloc");  // may be NULL on older glibc
}

// async-signal-safe integer to string (decimal)
static char *u64_to_str(uint64_t v, char *end) {
  if (v == 0) {
    *--end = '0';
    return end;
  }
  while (v) {
    *--end = '0' + (v % 10);
    v /= 10;
  }
  return end;
}

static void dump_stats(int signo) {
  (void)signo;
  char buf[256];
  char *p = buf;
  const char prefix[] = "HEAP current=";
  memcpy(p, prefix, sizeof(prefix) - 1);
  p += sizeof(prefix) - 1;
  char tmp[64];
  char *e = tmp + sizeof(tmp);
  uint64_t cur = atomic_load_explicit(&cur_bytes, memory_order_relaxed);
  uint64_t peak = atomic_load_explicit(&peak_bytes, memory_order_relaxed);
  char *s1 = u64_to_str(cur, e);
  memcpy(p, s1, (size_t)(e - s1));
  p += (e - s1);
  const char mid[] = " bytes, peak=";
  memcpy(p, mid, sizeof(mid) - 1);
  p += sizeof(mid) - 1;
  e = tmp + sizeof(tmp);
  char *s2 = u64_to_str(peak, e);
  memcpy(p, s2, (size_t)(e - s2));
  p += (e - s2);
  const char suffix[] = " bytes\n";
  memcpy(p, suffix, sizeof(suffix) - 1);
  p += sizeof(suffix) - 1;
  write(STDERR_FILENO, buf, (size_t)(p - buf));
}

__attribute__((constructor)) static void init_tracker(void) {
  resolve_syms();
  struct sigaction sa;
  memset(&sa, 0, sizeof(sa));
  sa.sa_handler = dump_stats;
  sigaction(SIGUSR1, &sa, NULL);
}

__attribute__((destructor)) static void fini_tracker(void) { dump_stats(0); }

static void account_add(size_t add) {
  if (!add) return;
  uint64_t cur =
      atomic_fetch_add_explicit(&cur_bytes, add, memory_order_relaxed) + add;
  uint64_t prev_peak = atomic_load_explicit(&peak_bytes, memory_order_relaxed);
  while (cur > prev_peak && !atomic_compare_exchange_weak_explicit(
                                &peak_bytes, &prev_peak, cur,
                                memory_order_relaxed, memory_order_relaxed)) {
    /* retry */
  }
}

static void account_sub(size_t sub) {
  if (!sub) return;
  atomic_fetch_sub_explicit(&cur_bytes, sub, memory_order_relaxed);
}

void *malloc(size_t size) {
  if (!real_malloc) resolve_syms();
  void *p = real_malloc(size);
  if (p) account_add(malloc_usable_size(p));
  return p;
}

void free(void *ptr) {
  if (!ptr) return;
  if (!real_free) resolve_syms();
  size_t sz = malloc_usable_size(ptr);
  account_sub(sz);
  real_free(ptr);
}

void *realloc(void *ptr, size_t size) {
  if (!real_realloc) resolve_syms();
  size_t old = ptr ? malloc_usable_size(ptr) : 0;
  void *np = real_realloc(ptr, size);
  if (np) {
    size_t neu = malloc_usable_size(np);
    if (neu >= old)
      account_add(neu - old);
    else
      account_sub(old - neu);
  } else {
    // realloc failure keeps old block intact
  }
  return np;
}

void *calloc(size_t nmemb, size_t size) {
  if (!real_calloc) resolve_syms();
  void *p = real_calloc(nmemb, size);
  if (p) account_add(malloc_usable_size(p));
  return p;
}

int posix_memalign(void **memptr, size_t alignment, size_t size) {
  if (!real_posix_memalign) resolve_syms();
  int r = real_posix_memalign(memptr, alignment, size);
  if (r == 0 && *memptr) account_add(malloc_usable_size(*memptr));
  return r;
}

void *aligned_alloc(size_t alignment, size_t size) {
  if (!real_aligned_alloc) resolve_syms();
  void *p = real_aligned_alloc ? real_aligned_alloc(alignment, size) : NULL;
  if (p) account_add(malloc_usable_size(p));
  return p;
}