/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef SRC_MULTITHREADDATA_H__
#define SRC_MULTITHREADDATA_H__

#include <err.h>
#include <pthread.h>
#include <semaphore.h>

#include <list>
#include <string>

#include "src/ReadPair.h"

template<class ITEM>
class ProtectedList {
 private:
  std::list <ITEM> items;
  pthread_mutex_t list_access;
  sem_t empty_slots;
  sem_t full_slots;
  size_t num_items;
  int slots;

 public:
  explicit ProtectedList(int _slots) :
  num_items(0),
    slots(_slots) {
      pthread_mutex_init(&list_access, NULL);
      sem_init(&empty_slots, 0, _slots);
      sem_init(&full_slots, 0, 0);
      num_items = 0;
    }

  void put(ITEM item) {
    sem_wait(&empty_slots);

    // Lock the list access
    pthread_mutex_lock(&list_access);

    // Add this record to the list
    items.push_back(item);

    // Release the lock
    pthread_mutex_unlock(&list_access);

    num_items++;

    // Flag consumer threads
    sem_post(&full_slots);
  }

  ITEM get() {
    // Wait for a semaphore to be avialble
    sem_wait(&full_slots);

    // Lock the list access
    pthread_mutex_lock(&list_access);

    // Get an element from the list, and remove it

    // Should never happen...
    if (items.empty())
      errx(1, "Internal error: Multithreaded " \
           "Protected list is empty (after getting a semaphore)");

    ITEM item = items.front();
    items.pop_front();
    num_items--;

    // Release the lock
    pthread_mutex_unlock(&list_access);
    sem_post(&empty_slots);
    return item;
  }

  void wait_for_all_slots() {
    for (int i = 0; i < slots; ++i) {
      sem_wait(&empty_slots);
    }
  }
};

class MultithreadData {
 private:
  ProtectedList<ReadPair*> items_to_process;
  ProtectedList<ReadPair*> items_to_output;

  int slots;
  pthread_mutex_t counter_mutex;
  size_t input_count;
  size_t output_count;

 public:
  explicit MultithreadData(int _slots);

  /* From the READER(Producer) to the Satellite
     processing threads (consumer) */
  void post_new_input_read(ReadPair* pRecord);

  /* Used in the Satellite processing threads */
  ReadPair* get_new_input();

  void wait_for_completed_input_processing();

  /* From the Satellite processing threads to the output-writer thread */
  void post_new_output_read(ReadPair *pRecord);

  /* Used in the Output-Writer thread */
  ReadPair* get_new_output();
  void wait_for_completed_output_processing();
  void increment_input_counter();
  void increment_output_counter();
  bool input_output_counters_equal();
};

#endif  // SRC_MULTITHREADDATA_H__
