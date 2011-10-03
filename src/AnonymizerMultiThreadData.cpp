/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include <err.h>
#include <iostream>

#include "AnonymizerMultiThreadData.h"

AnonymizerMultithreadData::AnonymizerMultithreadData(int _slots) :
  slots(_slots),
  items_to_process(_slots),
  input_count(0), output_count(0) {
  pthread_mutex_init(&counter_mutex,NULL);
}

void AnonymizerMultithreadData::increment_input_counter() {
  pthread_mutex_lock(&counter_mutex);
  input_count++;
  pthread_mutex_unlock(&counter_mutex);
}

void AnonymizerMultithreadData::increment_output_counter() {
  pthread_mutex_lock(&counter_mutex);
  output_count++;
  pthread_mutex_unlock(&counter_mutex);
}

bool AnonymizerMultithreadData::input_output_counters_equal() {
  bool equal;
  pthread_mutex_lock(&counter_mutex);
  std::cerr << "input_count = " << input_count << "  output_count = " << output_count << std::endl;
  equal =  ( output_count == input_count ) ;
  pthread_mutex_unlock(&counter_mutex);
  return equal;
}

void AnonymizerMultithreadData::post_new_input(const std::string& filename) {
  items_to_process.put(filename);
}

std::string AnonymizerMultithreadData::get_new_input() {
  return items_to_process.get();
}

void AnonymizerMultithreadData::wait_for_completed_input_processing() {
  items_to_process.wait_for_all_slots();
}

