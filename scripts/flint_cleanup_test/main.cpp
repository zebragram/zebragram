#include <iostream>
#include <omp.h>
#include <flint/fmpz.h>
#include <mutex>


// global lock
std::mutex print_lock;
uint64_t global_time = 0;
// RAII wrapper class
class FlintThreadManager {
public:
    FlintThreadManager() {
        // This constructor is called once per thread.
        // You could add initialization here if needed.
        // add lock for the print
        print_lock.lock();
        std::cout << "Thread " << omp_get_thread_num() 
                  << " has started, initializing flint." << std::endl;
        print_lock.unlock();
    }

    void manual_cleanup() {
        if (global_time == 0) return; // skip cleanup in the first phase
        print_lock.lock();
        std::cout << "Thread " << omp_get_thread_num() 
                  << " is manually calling flint_cleanup()." << std::endl;
        print_lock.unlock();
        flint_cleanup();
    }
    
    ~FlintThreadManager() {
        // This destructor is called only when the thread is destroyed.
        print_lock.lock();
        std::cout << "Thread " << omp_get_thread_num() 
                  << " is exiting, calling flint_cleanup()." << std::endl;
        print_lock.unlock();
        flint_cleanup();
    }
};

// Create one global, thread_local instance of the manager.
// C++ guarantees its destructor will be called when the thread exits.


void do_work() {
    // Thanks to the thread_local manager, we don't need manual cleanup here.
    thread_local FlintThreadManager flint_manager;
    fmpz_t x;
    fmpz_init_set_ui(x, omp_get_thread_num());
    fmpz_pow_ui(x, x, 10);
    // fmpz_clear is still needed for specific variables.
    fmpz_clear(x);
    // flint_manager.manual_cleanup(); // manually trigger cleanup at the end of work
}

int main() {
    global_time = 0;
    #pragma omp parallel
    {
        // The flint_manager is automatically created for each thread on first use.
        #pragma omp for
        for (int i = 0; i < 10; ++i) {
            do_work();
        }
    } // No cleanup needed here.
    global_time += 1;
    std::cout << "First phase done. Starting second phase." << std::endl;
    #pragma omp parallel for
    for (int i = 0; i < 10; ++i) {
        do_work();
    }

    std::cout << "Program finished. Cleanup will happen as threads exit." << std::endl;
    
    // The main thread also has its own flint_manager instance.
    // Its destructor will be called when main() returns.
    
    return 0;
}