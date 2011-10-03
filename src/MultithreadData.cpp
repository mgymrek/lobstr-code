/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include <err.h>
#include <iostream>

#include "MultithreadData.h"

MultithreadData::MultithreadData(int _slots) :
  slots(_slots),
  items_to_process(_slots), items_to_output(_slots),
  input_count(0), output_count(0) {
  pthread_mutex_init(&counter_mutex,NULL);
}

void MultithreadData::increment_input_counter() {
  pthread_mutex_lock(&counter_mutex);
  input_count++;
  pthread_mutex_unlock(&counter_mutex);
}

void MultithreadData::increment_output_counter() {
  pthread_mutex_lock(&counter_mutex);
  output_count++;
  pthread_mutex_unlock(&counter_mutex);
}

bool MultithreadData::input_output_counters_equal() {
  bool equal;
  pthread_mutex_lock(&counter_mutex);
  std::cerr << "input_count = " << input_count << "  output_count = " << output_count << std::endl;
  equal =  ( output_count == input_count ) ;
  pthread_mutex_unlock(&counter_mutex);
  return equal;
}

void MultithreadData::post_new_input_read(MSReadRecord* pRecord) {
  items_to_process.put(pRecord);
}

MSReadRecord* MultithreadData::get_new_input() {
  return items_to_process.get();
}

void MultithreadData::wait_for_completed_input_processing() {
  items_to_process.wait_for_all_slots();
}

void MultithreadData::post_new_output_read(MSReadRecord* pRecord) {
  items_to_output.put(pRecord);
}

MSReadRecord* MultithreadData::get_new_output() {
  return items_to_output.get();
}

void MultithreadData::wait_for_completed_output_processing() {
  items_to_output.wait_for_all_slots();
}
